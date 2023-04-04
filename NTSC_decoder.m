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

% find blanking pulse timing
pulseTiming = zeros(N,1);
pulseTiming(v > -.2) = 1;
pulseTiming = [diff(pulseTiming); 0];

% initialize loop vars
lineEnd = 0;
lineNo = 1;
frame = cell(1000,3);

% loop through every line in the whole signal
while (lineNo < 1000)
    % find start and end of line
    lineStart = min(t(pulseTiming == 1 & t > lineEnd));
    lineEnd = min(t(pulseTiming == -1 & t > lineStart));
    
    % break at end of recording
    if (isempty(lineStart) | isempty(lineEnd))
        break
    end

%     plot(t,[pulseTiming v], lineStart,0,'x', lineEnd,0,'x')
%     xlim([lineStart lineEnd])
%     ylim([-0.286 0.936])
    
    % grab line voltage
    lineV = v(t > lineStart & t < lineEnd)';
    
    % skip short lines
    if (length(lineV) < 800)
        continue
    end
    
    % trim line voltage
    lineV = lineV(1:7350);
    
    % store line data into frame
    frame(lineNo,:) = {lineStart lineEnd lineV};
    lineNo = lineNo + 1;
end

% display frame
imagesc(cell2mat(frame(:,3)))