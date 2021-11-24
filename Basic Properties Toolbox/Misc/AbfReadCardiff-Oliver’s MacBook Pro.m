function [myfilename] = AbfReadCardiff

%%  Author: O.G. Steele, 15.10.17
%%  Description
% AbfReadCardiff sequentially performs quiet version of abfload.m in order
% of .abf files in current directory
% In line with Cardiff Basic Properties SOP;
%   - file number 1 = GapFree protocol in Current Clamp mode (I=0)
%   - file number 2 = GapFree protocol in Current Clamp mode with the cell clamped to -70mV
%   - file number 3 = CurrentStep protocol in Current Clamp mode with sequentil injections of current from -10pA to +180pA in 10pA steps
%   - file number 4 = ActInact protocol in Voltage Clamp mode  with a holding potential of -70mV followed by 80ms steps from -120 to +80mV in increments of 5mV before being stepped for 200ms to the test potential of 0mV

%% Code
abf = dir('*.abf'); 
numfiles = length(abf);
mydata = cell(1, numfiles);

for k = 1:numfiles 
  mydata{k} = abfload(abf(k).name);
  myfilename = sprintf('file%d.txt', k);
end
GapFree = mydata{1,1};
GapFree70 = mydata{1,2};
CurrentStep = mydata{1,3};
ActInact = mydata{1,4};

end
 