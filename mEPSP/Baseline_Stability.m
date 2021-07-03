%% Baseline_Stab

% Author: O.G. Steele
% Last Update: 03.07.2021

% Description: A script used to take all of the events over a course of a
% recording and plot the average amplitude per wave

% To do: 
%   1. add in Vm compensation
%   2. save outputs
%   3. plot error bars and individual points

%% Load in event_data and baseline values

% load in event_data.phy with ephysIO
disp('Select event_data.phy')
[event_file, event_path] = uigetfile('*event_data.phy');
event_data = ephysIO(append(event_path,event_file));

% load in baseline values from txt files
idx = strfind(event_path,filesep);
base_path = event_path(1:idx(end-3));
cd(base_path)
base_file = dir('*baseline.txt');
baseline = table2array(readtable(base_file.name));


%% Split events into waves

% load in event_counts.txt
event_counts = table2array(readtable(fullfile(event_path,'event_counts.txt')));
cum_event_counts = cumsum(event_counts);
waves = size(event_counts,1);
n_events = sum(event_counts);
split_length = 10;

% display summary counts
disp(['Total number of waves = ', num2str(waves)])
disp(['Total number of events = ', num2str(n_events)])
disp(['Total average of frequency (Hz) = ', num2str(n_events/(waves*split_length))])

% bin events per minute, average and amplitude them
event_data_notime = event_data.array(:,2:end);
bins = 6:6:waves;
bins_1 = [1 bins];
for i = 1:size(bins,2)
    median_binned_events(:,i) = median(event_data_notime(:,bins_1(i):(cum_event_counts(bins(i)))),2);
    amp_binned_events(i) = max(median_binned_events(:,i));
end

% plot event amplitudes over wave number
figure; plot(bins/6,amp_binned_events*10^3); ylim([0,0.5]); xlabel('Time (minutes)'); ylabel('Amplitude (mV)');

