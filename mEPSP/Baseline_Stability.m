%% Baseline_Stab

% Author: O.G. Steele
% Last Update: 03.07.2021

% Description: A script used to take all of the events over a course of a
% recording and plot the average amplitude per wave

% To do: 
%   1. add in Vm compensation

%% Load in event_data and baseline values

% clean environment
clear; close all

% load in event_data.phy with ephysIO
disp('Select event_data.phy')
[event_file, event_path] = uigetfile('*event_data.phy');
event_data = ephysIO(append(event_path,event_file));
cd(event_path)

% load in baseline values from txt files
filenameVar = char(table2array(readtable('..\..\filepaths.txt')));
basename = [filenameVar(3:end-5) '_baseline.txt'];
cd ..\..\..
baseline = table2array(readtable(basename));


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
event_amps = max(event_data_notime(:,:));
for i = 1:size(bins,2)
    median_binned_events(:,i) = median(event_data_notime(:,bins_1(i):(cum_event_counts(bins(i)))),2);
    %mean_binned_events(:,i) = mean(event_data_notime(:,bins_1(i):(cum_event_counts(bins(i)))),2);
    amp_median_binned_events(i) = max(median_binned_events(:,i));
    %amp_mean_binned_events(i) = max(mean_binned_events(:,i));
    amp_stderror(i) = std( event_amps(:,bins_1(i):(cum_event_counts(bins(i)))) ) / sqrt( length( event_amps(:,bins_1(i):(cum_event_counts(bins(i)))) ))*10^3;
    binned_baseline(i) = median(baseline(bins_1(i):bins(i),2));
    base_stderror(i) = std( baseline(bins_1(i):bins(i),2) / sqrt( length ( baseline(bins_1(i):bins(i),2))));   
end

time = 1:size(bins,2);

% plot holding potential
figure; plot(baseline(:,1)/60,baseline(:,2),'LineWidth',2); xlabel('Time (minutes)'); ylabel('Membrane potential (mV)');
box off; set(gcf,'color','w'); set(gca,'LineWidth',2); title('Membrane potential over time');

% plot event amplitudes (mV) over wave number
figure; plot(time,amp_median_binned_events*10^3,'-o','MarkerFaceColor','b'); ylim([0,0.5]); xlabel('Time (minutes)'); ylabel('Amplitude (mV)');
hold on; errorbar(time,amp_median_binned_events*10^3,amp_stderror,'Color','black','LineStyle','none'); box off; set(gcf,'color','w'); set(gca,'LineWidth',2);
title('Median event amplitude +/- SEM');
%figure; plot(bins/6,amp_mean_binned_events*10^3); ylim([0,0.5]); xlabel('Time (minutes)'); ylabel('Amplitude (mV)');

saveas(gcf,'baseline_stability.pdf')
