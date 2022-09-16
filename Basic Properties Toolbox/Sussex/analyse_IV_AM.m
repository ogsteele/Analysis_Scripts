[file,path] = uigetfile('*.ma'); % get the file of interest
clampfile = fullfile(path,file); % generate the full file path of the file
SI = ephysIO({clampfile,2}); % load in the current data
x = SI.array(:,1); % x is time here, but not in seconds, data points
Y = SI.array(:,2:end); % y is the array of current data
base = mean (mean (Y(4001:16001,:))); % determine the baseline (intra-step)
N = size(SI.array,2) - 1; % get the number of waves in the array (minus time)
lo = dsearchn(x,0.5) + 1; % start of the test pulse
hi = dsearchn(x,1.5) - 1; % end of the test pulse
steps_all = zeros (1,N); % preallocate a blank series of steps
for i = 1:N
   steps_all(i) = mean (Y(lo:hi,i)); % fill the steps with the mean of each step
end
steps = fix(([steps_all(1) steps_all(end)] - base) * 1e+12) % round to zero, baseline subtract and put into pA
SV = ephysIO({clampfile,1}); % load the voltage data
output = IStep(SV,clampfile,steps); % run the Istep function

% note to OGS: look to integrate the above into Istep