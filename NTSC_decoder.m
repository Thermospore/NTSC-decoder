% using a plain old .m so git can see changes
clear all
close all

bonusLines = 0; % note: semi-temp fix. accounts for how n64 doesn't have half width vsync pulses...
rfCH = -1; % used to set rf demod (-1 = no demod)

%% set parameters

sat = 1.3; % saturation
bri = 1.3; % brightness
outputScale = 3; % gif output res is multiplied by this factor
speedScale = .5; % gif framerate is multiplied by this factor

% uncomment a recording to view it

% note: figure out why the ps1/2 recordings are so noisy...
%       not sure if it's specific to the 10frame preset or the ps2
% note: add a pure noise mode :)
% note: implement interlaced video

% recordingName = "2frame_attract"; dcOffset = -0.1329; lowPulse = -0.7525; bonusLines = 6;
% recordingName = "2frame_mystery"; dcOffset = -0.1918; lowPulse = -0.8077; bonusLines = 6;
% recordingName = "2frame_n64logo"; dcOffset = -0.1796; lowPulse = -0.7944; bonusLines = 6;
% recordingName = "manyframe_mystery"; dcOffset = -0.1896; lowPulse = -0.8068; bonusLines = 6;
recordingName = "10frame_mystery"; dcOffset = -0.0826; lowPulse = -0.6682;
% recordingName = "10frame_ps1startup"; dcOffset = 0.0308; lowPulse = -0.5623;
% recordingName = "10frame_ps2startup"; dcOffset = 0.0666; lowPulse = -0.5283;
% recordingName = "2frame_f0menu"; rfCH = 3; % this one is non functional
% recordingName = "1frame_f0jump"; dcOffset = 1.0950e-3; lowPulse = 1.3347e-3; rfCH = 3;
% recordingName = "10frame_c2intro"; dcOffset = 0.0131; lowPulse = -0.5793;
% recordingName = "10frame_c2intro_take2"; dcOffset = 0.0025; lowPulse = -0.5917;
% recordingName = "10frame_c2menu"; dcOffset = -0.1612; lowPulse = -0.7528;
% recordingName = "10frame_75load_c2menu"; dcOffset = -0.2764; lowPulse = -0.5712;
% recordingName = "manyframe_75load_mystery"; dcOffset = -0.0515; lowPulse = -0.3446;

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
zTotal_line_period = 63.556e-6; % (Mono: 63.492ƒÊs)
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
linesPerFrame = 254+bonusLines; % (# full length scanlines in 240p)
f_SC = (315/88)*1e6; % Color subcarrier freq

%% import data
disp("importing data")

load("scope recordings\" + recordingName + ".mat")

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

% lazy/hardcoded rf demod
if (rfCH ~= -1)
    VV = fftshift(fft(v));
    VV = fftshift([VV(5225047:end); VV(1:5225047-1)]);
    VV(abs(ff)> 4.4e6) = 0;
    v = abs(ifft(ifftshift(VV)));
end

% fix scaling
% note: maybe implement autoscaling... can use that vsync pulse finder code I wrote
v = (v-dcOffset)/(lowPulse-dcOffset)*zSync_tip_level;

%% separate signals
disp("separating signals")

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
% note: maybe normalize to expected amplitude of color burst?

% extract timing info
% note: consider applying LPF to v first, to prevent extra pulses due to noise?
%       but I think they get ignored anyway, when we toss out short lines?
%       wait that would throw off vsync pulse detection tho, if an extra
%       pulse was mixed in there...
pulv = zeros(N,1);
if (rfCH == -1)
    pulv(filtfilt(h,1,v) > zSync_tip_level/2) = 1;
else
    % extra noise filtering for rf demod case
    pulv(13*filtfilt(h,1,filtfilt(h,1,v)) > zSync_tip_level/2) = 1;
end
pulv = [diff(pulv); 0];

lineStarts = t(pulv == 1);
lineEnds = t(pulv == -1);
if (lineEnds(1) < lineStarts(1))
    lineEnds = lineEnds(2:end);
end
lineSegs = [lineStarts(1:length(lineEnds)) lineEnds];

%% main render loop
disp("running main render loop")

% loop through every line in the whole signal
lineNo = 1;
lineNo_firstFullFrame = -1;
lineEnd = 0;
lineLengthHist = [nan nan nan];
lineN = floor((zTotal_line_period-zLine_sync_pulse_width - .05*zFront_porch)/T);
longFrame = zeros(length(lineSegs),lineN,3); % pre-allocate (guess length)
for i = 1:length(lineSegs)
    
    % get start and end of line
    lineStart = lineSegs(i,1);
    lineEnd = lineSegs(i,2);
    
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
    index = ceil(lineStart/T):ceil(lineStart/T)+lineN-1;
    linelum = lumv(index)';
    linechr = chrc(index).'; % DON'T FORGET DOT!!!
    
    % note: fix line voltage dc coupling?
    %       and remove sync pulse here instead? (ie slice off <0 stuff)
    %       US vs JP black level?
    
    % set color burst to 180deg
    index = ceil((zTime_reference_point_to_burst_start-zLine_sync_pulse_width + 1/6*zSubcarrier_burst_duration)/T):...
            floor((zTime_reference_point_to_burst_start-zLine_sync_pulse_width + 5/6*zSubcarrier_burst_duration)/T);
    cbPhase = mean(unwrap(angle(linechr(index))));
    linechr = linechr*exp(j*(pi-cbPhase));
    
    % convert to HSV
    hsv_h = angle(linechr*exp(j*deg2rad(-110))); % HSV and NTSC have 110deg hue offset
    hsv_h(hsv_h<0) = hsv_h(hsv_h<0)+2*pi;
    hsv_h = hsv_h/2/pi;
    
    hsv_s = abs(linechr)/(zMaximum_excursion_with_chroma-zPeak_white_level);
    hsv_s = hsv_s*sat;
    hsv_s = min(1,hsv_s);
    
    hsv_v = linelum/zPeak_white_level;
    hsv_v = hsv_v*bri;
    hsv_v = min(1,hsv_v);
    
    % store line data into big long frame
    longFrame(lineNo, :, :) = cat(3,hsv_h,hsv_s,hsv_v);
    
    % update lineNo
    lineNo = lineNo + 1;
end

%% final outputs
disp("slicing final outputs")

% trim excess from pre-allocation
longFrame = longFrame(1:lineNo-1,:,:);

% convert from HSV to RGB
longFrame = hsv2rgb(longFrame);

% display big long frame
imagesc(longFrame)

% slice into complete frames for gif
for i = 1:floor((size(longFrame,1)-lineNo_firstFullFrame)/linesPerFrame)
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