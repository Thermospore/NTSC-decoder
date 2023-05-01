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

%     plot(t,[pulv v], lineStart,0,'x', lineEnd,0,'x')
%     xlim([lineStart lineEnd])
%     ylim([-0.286 0.936])
    
    % grab line voltage
    linev = v(t > lineStart & t < lineEnd)';
    
    % skip short lines
    if (length(linev) < 800)
        continue
    end
    
    % trim line voltage
    linev = linev(1:7350);
    
    % store line data into frame
    frame(lineNo,:) = {lineStart lineEnd linev};
    lineNo = lineNo + 1
end

% display frame
imagesc(cell2mat(frame(:,3)))