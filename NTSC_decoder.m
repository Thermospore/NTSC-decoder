% using a plain old .m so git can see changes
clear all
close all

%% set parameters

recordingName = "2frame_attract"; dcOffset = -0.1330;lowPulse = -0.7527 - dcOffset;
% recordingName = "2frame_mystery";
% recordingName = "2frame_n64logo";
% recordingName = "manyframe_mystery";

%% import data

load("scope recordings\" + recordingName + ".mat")

% fix scaling
v = (v-dcOffset)/lowPulse*(-0.286);

% calculate signal specs
N = length(t);
T = mean(diff(t));

%% do stuff

% extract luminance and chrominance info
VV = fftshift(fft(v));
ff = linspace(-1/T/2,1/T/2,N)';

combPeriod = 15.734265734265e3;
combHH = -cos(2*pi/combPeriod*ff) /2 + .5;

bandHH = 2*tripuls(ff,2*2e6);
bandHH = -min(bandHH,1) + 1;

chrHH = combHH.*bandHH;
lumHH = 1 - chrHH;

chrVV = chrHH.*VV;
lumVV = lumHH.*VV;

chrv = real(ifft(ifftshift(chrVV)));
lumv = real(ifft(ifftshift(lumVV)));

% demodulate chrominance
fsc = 3579545.45454545;

L = 100;
f = [0 fsc*.6 fsc*.8 fsc fsc*1.2 1/(2*T)]*2*T;
m = [1 1 0 0 0 0];
h = fir2(L,f,m);

chrr = chrv.*sin(2*pi*fsc*t);
chrr = filtfilt(h,1,chrr);
chrj = chrv.*cos(2*pi*fsc*t);
chrj = filtfilt(h,1,chrj);

chrc = chrr + j*chrj;

% extract timing info
pulv = zeros(N,1);
pulv(v > -.2) = 1;
pulv = [diff(pulv); 0];

% initialize loop vars
lineEnd = 0;
lineNo = 1;
frame = cell(1000,3);

% loop through every line in the whole signal
while (lineNo < 1000)
    % find start and end of line
    lineStart = min(t(pulv == 1 & t > lineEnd));
    lineEnd = min(t(pulv == -1 & t > lineStart));
    
    % break at end of recording
    if (isempty(lineStart) | isempty(lineEnd))
        break
    end

%     if(lineNo == 170)
%         plot(t,[v pulv lumv abs(chrc) angle(chrc)/pi], lineStart,0,'x', lineEnd,0,'x')
%         xlim([lineStart lineEnd])
%         ylim([-0.286 0.936])
%         break
%     end
    
    % grab line voltage
    linelum = lumv(t > lineStart & t < lineEnd)';
    linechr = chrc(t > lineStart & t < lineEnd).'; % DON'T FORGET DOT!!!
    
    % skip short lines
    if (length(linelum) < 800)
        continue
    end
    
    % set color burst to 180deg
    index = (t > lineStart + (5.3-4.7 + 2.5*1/6)*1e-6) &...
            (t < lineStart + (5.3-4.7 + 2.5*5/6)*1e-6);
    cbPhase = mean(angle(chrc(index)));
    linechr = linechr*exp(j*(pi-cbPhase));
    
    % extract U/V from chrominance
    saturation = 1.5;
    lineU = saturation*real(linechr);
    lineV = saturation*imag(linechr);
    
    % decode to RGB
    clear RGB
    RGB(:,:,1) = linelum + 1.13*lineU;
    RGB(:,:,2) = linelum + -0.575*lineU - 3.94*lineV;
    RGB(:,:,3) = linelum + 2.03*lineV;
    
    % make line length uniform
    RGB = RGB(:, 1:7350, :);
    
    % store line data into frame
    frame(lineNo,:) = {lineStart lineEnd RGB};
    lineNo = lineNo + 1
end

% display frame
gain = 1.2;
imagesc(gain*cell2mat(frame(:,3)))