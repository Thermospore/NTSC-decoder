% using a plain old .m so git can see changes
clear all
close all

%% set parameters

sat = 1.3; % saturation
bri = 1.3; % brightness
outputScale = 3; % gif output res is multiplied by this factor
speedScale = .15; % gif framerate is multiplied by this factor

% need to implement autoscaling...
% can use that vsync pulse finder code I wrote
recordingName = "2frame_attract"; dcOffset = -0.1330;lowPulse = -0.7527 - dcOffset;
% recordingName = "2frame_mystery";
% recordingName = "2frame_n64logo";
% recordingName = "manyframe_mystery";

%% define NTSC constants
% source: https://web.archive.org/web/20170614080536/http://www.radios-tv.co.uk/Pembers/World-TV-Standards/Line-Standards.html

% 525 Lines
% everything is wrt 50% amplitude points
% rise times are 10% to 90%

% the z is just so they don't clutter up workspace
% (probably a sleeker way to do this...)

% Amplitudes (volts)
zMaximum_excursion_with_chroma = 936e-3;
zPeak_white_level = 714e-3;
zPositive_peaks_of_burst = 143e-3;
zBlack_level = 0; % 54 mV outside Japan
zBlanking_level = 0;
zNegative_peaks_of_burst = -143e-3;
zMinimum_excursion_with_chroma = -164e-3;
zSync_tip_level = -286e-3;

% Horizontal Timing (seconds unless specified)
zLine_frequency = (2.25/143)*1e6; % Hz (Mono: 15.750kHz)
zTotal_line_period = 63.556e-6; % (Mono: 63.492��s)
zActive_line_period = 52.9e-6;
zLine_blanking_period = 10.7e-6;
zLine_blanking_rise_time = 140e-9;
zFront_porch = 1.5e-6;
zLine_sync_pulse_width = 4.7e-6;
zLine_sync_rise_time = 140e-9;
zTime_reference_point_to_burst_start = 5.3e-6; % (19 cycles)
zSubcarrier_burst_phase_wrt_time_reference_point = 0; % degrees
zSubcarrier_burst_duration = 2.5e-6; % (9 cycles)
zBurst_envelope_rise_time = 300e-9;
zTime_reference_point_to_line_blanking_end = 9.2e-6;

% Vertical Timing
zField_frequency = 60/1.001; % Hz (Mono: 60Hz)
zTotal_field_period = 16.6833e-3; % sec (Mono: 16.667ms)
zTotal_frame_period = 33.3667e-3; % sec (Mono: 33.333ms)
zActive_field_period = 242.5; % lines
zField_blanking_period = 20; % lines (N lines plus line blanking)
zPre_equalising_duration = 3; % lines (6 narrow pulses)
zPre_equalising_pulse_width = 2.30e-6; % sec
zEqualising_pulse_rise_time = 140e-9; % sec
zField_sync_duration = 3; % lines (6 broad pulses)
zField_sync_broad_pulse_width = 27.1e-6; % sec
zField_serration_pulse_width = 4.70e-6; % sec
zField_pulse_rise_time = 140e-9; % sec
zPost_equalising_duration = 3; % lines (6 narrow pulses)
zPost_equalising_pulse_width = 2.30e-6; % sec

% misc
linesPerFrame = 260; % (# full length scanlines in 240p)
f_SC = (315/88)*1e6; % Color subcarrier

%% import data

load("scope recordings\" + recordingName + ".mat")

% fix scaling
v = (v-dcOffset)/lowPulse*zSync_tip_level;

% calculate time/freq vectors
N = length(t);
T = mean(diff(t));
NT = N*T;

F = 1/NT;
NF = 1/T;

f = transpose(0:F:NF-F);

tt = fftshift(t);
tt(tt > tt(end)) = tt(tt > tt(end)) - NT;
ff = fftshift(f);
ff(ff > ff(end)) = ff(ff > ff(end)) - NF;

%% do stuff

% extract luminance and chrominance info
VV = fftshift(fft(v));

combHH = -cos(2*pi/zLine_frequency*ff) /2 + .5;

bandHH = 2*tripuls(ff,2*2e6);
bandHH = -min(bandHH,1) + 1;

chrHH = combHH.*bandHH;
lumHH = 1 - chrHH;

chrVV = chrHH.*VV;
lumVV = lumHH.*VV;

chrv = real(ifft(ifftshift(chrVV)));
lumv = real(ifft(ifftshift(lumVV)));
lumv = max(0,lumv); % remove sync pulses

% demodulate chrominance
L = 100;
f = [0 f_SC*.6 f_SC*.8 f_SC f_SC*1.2 1/(2*T)]*2*T;
m = [1 1 0 0 0 0];
h = fir2(L,f,m);

chrr = chrv.*sin(2*pi*f_SC*t);
chrr = filtfilt(h,1,chrr);
chrj = chrv.*cos(2*pi*f_SC*t);
chrj = filtfilt(h,1,chrj);

chrc = chrr + j*chrj;
chrc = chrc/max(abs(chrc))*max(abs(chrv)); % (lazily and sloppily) restore amplitude scale

% extract timing info
% consider applying LPF to v first, to prevent extra pulses due to noise?
% but I think they get ignored anyway, when we toss out short lines
pulv = zeros(N,1);
pulv(v > (zMinimum_excursion_with_chroma + zSync_tip_level)/2) = 1;
pulv = [diff(pulv); 0];

% loop through every line in the whole signal
% go back and pull out things that only need to be calc'd once
lineNo = 1;
lineNo_firstFullFrame = -1;
lineEnd = 0;
lineLengthHist = [nan nan nan];
while (1 < 2)
    % find start and end of line
    lineStart = min(t(pulv == 1 & t > lineEnd));
    lineEnd = min(t(pulv == -1 & t > lineStart));
    
    % break at end of recording
    if (isempty(lineStart) | isempty(lineEnd))
        break
    end

%     if(lineNo == 164)
%         plot(t,[v pulv lumv abs(chrc) angle(chrc)/pi], lineStart,0,'x', lineEnd,0,'x')
%         xlim([lineStart lineEnd])
%         ylim([-0.286 0.936])
%         break
%     end
    
    % check for start of first full frame, if we haven't found it already
    lineLength = lineEnd - lineStart;
    if (lineNo_firstFullFrame == -1)
        
        % update history of last 3 line lengths
        lineLengthHist = [lineLengthHist(2:3) lineLength];

        % if they are all basically equal to zField_serration_pulse_width,
        % then we are at the start of the frame
        if (lineLengthHist > .95*zField_serration_pulse_width &...
            lineLengthHist < 1.05*zField_serration_pulse_width)
        
            lineNo_firstFullFrame = lineNo;
        end
    end
    
    % skip non-image lines
    if (lineLength < .95*(zTotal_line_period - zLine_sync_pulse_width))
        continue
    end
    
    % grab line voltage
    linelum = lumv(t > lineStart & t < lineEnd)';
    linechr = chrc(t > lineStart & t < lineEnd).'; % DON'T FORGET DOT!!!
    
    % set color burst to 180deg
    index = (t > lineStart + zTime_reference_point_to_burst_start-zLine_sync_pulse_width + 1/6*zSubcarrier_burst_duration) &...
            (t < lineStart + zTime_reference_point_to_burst_start-zLine_sync_pulse_width + 5/6*zSubcarrier_burst_duration);
    cbPhase = mean(angle(chrc(index)));
    linechr = linechr*exp(j*(pi-cbPhase));
    
    % decode to RGB
    hsv_h = angle(linechr*exp(j*deg2rad(-110))); % HSV and NTSC have 110deg hue offset
    hsv_h(hsv_h<0) = hsv_h(hsv_h<0)+2*pi;
    hsv_h = hsv_h/2/pi;
    
    hsv_s = abs(linechr)/(zMaximum_excursion_with_chroma-zPeak_white_level);
    hsv_s = hsv_s*sat;
    hsv_s = min(1,hsv_s);
    
    hsv_v = linelum/zPeak_white_level;
    hsv_v = hsv_v*bri;
    hsv_v = min(1,hsv_v);
    
    RGB = hsv2rgb(hsv_h,hsv_s,hsv_v);
    
    % make line length uniform
    RGB = RGB(:, 1:floor((zTotal_line_period-zLine_sync_pulse_width - .05*zFront_porch)/T), :);
    
    % store line data into big long frame
    longFrame(lineNo, :, :) = RGB;
    
    % update lineNo and display progress
    lineNo = lineNo + 1;
    totalProgress = [lineNo lineEnd/max(t)]
end

% display big long frame
imagesc(longFrame)

% slice into complete frames for gif
for i = 1:floor(size(longFrame,1)/linesPerFrame)
    frame(:,:,:,i) = longFrame(lineNo_firstFullFrame + (i-1)*linesPerFrame + (0:linesPerFrame-1),:,:);
end

% gif output
filename = "output.gif";
for idx = 1:size(frame,4)
    outputRes = [linesPerFrame*outputScale round(linesPerFrame/3*4*outputScale)]; % somewhat arbitrary aspect...
    outputFrame = imresize(frame(:,:,:,idx),outputRes,"nearest");
    [A,map] = rgb2ind(outputFrame,256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",1/zField_frequency/speedScale);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",1/zField_frequency/speedScale);
    end
end

winopen output.gif