function [output,Penn_SlopeVal,Max_Amp_Val] = fEPSP_Base
% Author: O.G. Steele
% Date of Creation: 21.02.19
% Updated: 19.08.21

% Change notes: ephysIO update has changed functionality, switching to
% h5read (19.08.21)
%% Arguments
% LPF = Low pass frequency value in Hz

%% Notes
% Make ACQ4_Base the master, with input arguements to select type
% Allow the script to loop
% reorganise final structure
% PeakAmp?

%% Analysis
% Clear windows
close all

% Change starting folder to MATLAB/PhD
%if ispc
    %start = cellstr('C:\Users\');
    %middle = cellstr(getenv('username'));
    %finish = cellstr('\OneDrive - University of Sussex\MATLAB\PhD');
    %filepath = fullfile(start,middle,finish);
    %cd(filepath{1});
%elseif ismac
    %start = cellstr('/Users/');
    %splitpwd = split(pwd,'/');
    %user = cellstr(splitpwd(3));
    %folder = cellstr('/OneDrive - University of Sussex/MATLAB/PhD');
    %macpath = fullfile(start,user,folder);
    %cd(char(macpath));
%end

%% ephysIO file selection
% Time_Diff and meta taken from first recording in first file path

path = uigetdir;
S = ephysIO({[path,'/000/Clamp2.ma'],1}); % loads channel 1 by default
Meta = readMeta([path,'/000/Clamp2.ma']); % reads metadata on the first clamp2.ma file, should be consistent across all

Time = S.array(:,1);
Time_Diff = S.xdiff;

if path == 0
	return;
end
topLevelFolder = path;

% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);

% Parse into a cell array.
listOfFolderNames = split(allSubFolders,';'); % single use of split over strtok
listOfFolderNames = listOfFolderNames(2:end,:); % remove the first as not a subfolder
listOfFolderNames = listOfFolderNames(~cellfun('isempty',listOfFolderNames)); % remove any empty values
numberOfFolders = length(listOfFolderNames); % count the number of folders


% Change directory to first sweep 
% Note, not first folder as that's the master folder
for i = 2:numberOfFolders
    cd(char(listOfFolderNames(i,:)));
    Data = h5read('Clamp2.ma','/data');
    Raw_Data(:,i-1) = Data(:,1);
end
%%
num_Sweeps = size(Raw_Data,2);

% adjust sweeps to zero and smooth to compensate for noise
basemean = mean(Raw_Data(1:3500,:));
zeroed_trace = (Raw_Data(:,:) - basemean) * 1000; % and unit converted to mV

% Apply a low pass 1000 Hz binomial filter to the data
YF = filter1(zeroed_trace, Time, 0, 1000);

% plot the sweeps and get points pre and post field potential
plot(zeroed_trace(:,:));
title('Select fEPSP Window')
ylabel('Amplitude (mV)')
xlabel('Number of data points')
%ylim([-1 1.5]);
%xlim([0 2000]);

% allow user to zoom and then hit enter when at region to exit zoom
zoom on
waitfor(gcf, 'CurrentCharacter', char(13))
zoom reset
zoom off

hold on
[x,~] = ginput;
close

% Seperate out window of interest
fEPSP_start = round(x(1)+3500);
fEPSP_end = round(x(2)+3500);
fEPSP_Window = YF((fEPSP_start:fEPSP_end),:);
num_points = size(fEPSP_Window,1);
fEPSP_Window_Time = linspace(0,Time_Diff * num_points,num_points);

% Index of maximum amplitude, and incidentally the peak amplitude of each
% stimulation intensity, preallocated for speed
Max_Amp_Val = zeros(1,num_Sweeps);
Max_Amp_Ind = zeros(1,num_Sweeps);

for o = 1:num_Sweeps
    [Max_Amp_Val(o),Max_Amp_Ind(o)] = min(fEPSP_Window(:,o));
end

%% Randall Slope

SlopeFilter = filter1(zeroed_trace(fEPSP_start:fEPSP_end,:), fEPSP_Window_Time, 0, 250);

% index of 30% and 70% of maximum amplitudes
I_slope_30 = round(Max_Amp_Ind*0.3);
I_slope_70 = round(Max_Amp_Ind*0.7);

% calculate fEPSP slope between 30% and 70% peak amplitude
Randall_slope = zeros(1,num_Sweeps);
for k = 1:num_Sweeps
    Randall_slope(k) = (mean(gradient(SlopeFilter((I_slope_30:I_slope_70),k))));
end

%% Penn Slope
% Calculates slope as the peak of the first derivative of the slope.
% Slope determined as start to maximum amplitude
% Aim to make the peak of the first derivative between 0.3 and 0.7 of the
% slope - should reduce noise! Merging Penn and Randall slope detection.

% Calculates slope in mV/S

SlopeFilter = filter1(zeroed_trace(fEPSP_start:fEPSP_end,:), fEPSP_Window_Time, 0, 250);

for x = 1:num_Sweeps
    SlopeFilter(:,x) = SlopeFilter(:,x) - mean(SlopeFilter(1:5,x));
end

Penn_SlopeVal = zeros(1,num_Sweeps);
Max = zeros(1,num_Sweeps);
ind = zeros(1,num_Sweeps);
SlopeVal_Index = zeros(1,num_Sweeps);

for x = 1:num_Sweeps
    [Max(x),ind(x)] = min(SlopeFilter(:,x));
end

I_Slope_30 = ceil(ind*0.2);
I_Slope_70 = ceil(ind*0.8)+2; % +2 to account for no negative slope error

for x = 1:num_Sweeps
    [dydx, ~, ~] = ndiff(SlopeFilter(I_Slope_30(x):I_Slope_70(x),x), fEPSP_Window_Time(I_Slope_30(x):I_Slope_70(x)));
    [Penn_SlopeVal(x), SlopeVal_Index(x)] = min(dydx);
end


%% Plotting

% Plot fEPSP Traces
figure('units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(2,3,1)
for l = 1:num_Sweeps
    %base = mean(fEPSP_Window(1:5,l));
    plot(fEPSP_Window_Time * 1000,(fEPSP_Window(:,l)))
    hold on
    for i = 1:num_Sweeps
        plot(fEPSP_Window_Time(Max_Amp_Ind(i))*1000,Max_Amp_Val(i),'o','Color','r','LineWidth',3) 
    end
end
ylabel('Amplitude (mV)','fontsize',16);
xlabel('Time (ms)','fontsize',16);
title('1 kHz LPF fEPSP + Max. Amplitude','fontsize',18)
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
hold off

% Plot 30-70 overlay
subplot(2,3,2)
plot(fEPSP_Window_Time * 1000,SlopeFilter,'r');
hold on
for x = 1:num_Sweeps
    hold on
    plot(fEPSP_Window_Time(I_Slope_30(x):I_Slope_70(x)) * 1000, SlopeFilter(I_Slope_30(x):I_Slope_70(x),x),'b');
    plot(fEPSP_Window_Time(I_Slope_30(x) + SlopeVal_Index(x)) * 1000, SlopeFilter((I_Slope_30(x) + SlopeVal_Index(x)),x),'o','Color','k','LineWidth',3);
end
box off
xlabel('Time (ms)')
ylabel('Amplitude (mV)')
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
title('250 Hz LPF fEPSP, 20-80% Max Amp, Max dydx')
hold off

% Plot Median Trace
% Notes: Plots median of the first 6th, middle 6th and final 6th
% In the case of a 60 minute recording, 1/6th = 10 minutes = 20 sweeps

subplot(2,3,3)
    median_trace_early = (median(fEPSP_Window(:,1:20),2));
    base_median_early = mean(median_trace_early(1:5));
    zeroed_median_early = (median_trace_early - base_median_early);
        plot(fEPSP_Window_Time * 1000 ,zeroed_median_early,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineWidth',4,'DisplayName','Early Median')
        hold on
        array = fEPSP_Window(:,1:20)*-1;
        for i = 1:size(array,1)
             stdev(i) = std(array(i,:));
        end
        %errorbar(fEPSP_Window_Time * 1000 ,zeroed_median_early,stdev,'Color',[0 0.4470 0.7410])
        shadedErrorBar(fEPSP_Window_Time * 1000 ,zeroed_median_early,stdev,'lineprops','b')
    median_trace_late = (median(fEPSP_Window(:,end - 20:end),2));
    base_median_late = mean(median_trace_late(1:5));
    zeroed_median_late = (median_trace_late - base_median_late);
plot(fEPSP_Window_Time * 1000 ,zeroed_median_late,'LineWidth',4,'Color',[0.8500, 0.3250, 0.0980],'DisplayName','Late Median')
     array = fEPSP_Window(:,end - 20:end)*-1;
            for i = 1:size(array,1)
                 stdev(i) = std(array(i,:));
            end
            %errorbar(fEPSP_Window_Time * 1000 ,zeroed_median_late,stdev,'Color',[0.8500, 0.3250, 0.0980])
            shadedErrorBar(fEPSP_Window_Time * 1000 ,zeroed_median_late,stdev,'lineprops','r')
legend('Early Median','Early Median SEM','Late Median','Late Median SEM','Location','Northeast','linewidth',1)
title('Median fEPSP')
box off
xlabel('Time (ms)')
ylabel('Median fEPSP Traces')
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)

% Plot amplitude
subplot(2,3,4)
plot(linspace(0,num_Sweeps/2,num_Sweeps),(Max_Amp_Val * -1));
hold on
Max_Amp_Val_movmean30 = movmean(Max_Amp_Val * -1,num_Sweeps/5);
plot(linspace(0,num_Sweeps/2,num_Sweeps),Max_Amp_Val_movmean30)
ylabel('Amplitude (mV)','fontsize',16);
xlabel('Time (minutes)','fontsize',16);
title('Maximum Amplitude (mV)','fontsize',18);
ylim([0 (max(Max_Amp_Val * -1) * 1.5)])
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
set(gcf, 'Position',  [100, 100, 1500, 400])
legend('Amplitude (mV)','Moving Average','Location','Northwest','linewidth',1)

hold off

% Plot Penn slope
subplot(2,3,5)
plot(linspace(0,num_Sweeps/2,num_Sweeps),Penn_SlopeVal * -1);
hold on
SlopeVal_movmean20 = movmean(Penn_SlopeVal * -1,num_Sweeps/5);
plot(linspace(0,num_Sweeps/2,num_Sweeps),SlopeVal_movmean20)
ylabel('fEPSP slope (mV/s)','fontsize',16);
xlabel('Time (minutes)','fontsize',16);
title('Penn fEPSP slope (mV/s)','fontsize',18);
ylim([0 (max(Penn_SlopeVal * -1) * 1.5)])
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
legend('Slope (mV/s)','Moving Average','Location','Northwest','linewidth',1)
hold off

% Plot Randall slope
subplot(2,3,6)
plot(linspace(0,num_Sweeps/2,num_Sweeps),(Randall_slope * -100000)/2.5)
hold on
Slope_movmean20 = movmean((Randall_slope * -100000)/2.5,num_Sweeps/5);
plot(linspace(0,num_Sweeps/2,num_Sweeps),Slope_movmean20)
box off
title('Randall fEPSP Slope(mV/s)')
xlabel('Time (minutes)')
ylabel('Slope Value (mV/s)')
ylim([0 (max((Randall_slope * -100000)/2.5) * 1.5)])
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
legend('Slope (mV/s)','Moving Average','Location','Northwest','linewidth',1)

%% Create output structure
%output.HalfMax = HM;
output.ephysIO = S;
output.Meta = Meta;
output.Data.Max_Amp_Val = Max_Amp_Val';
output.Data.Max_Amp_Ind = Max_Amp_Ind;
output.Data.Slope = Randall_slope;
output.Data.Raw_Data = zeroed_trace;
output.Data.filtered_data = YF;
output.Data.fEPSP_Window = fEPSP_Window;
output.Data.fEPSP_Window_Time = fEPSP_Window_Time;
output.Data.SlopeVals = Penn_SlopeVal';
output.Data.SlopeVal_Ind = SlopeVal_Index;
output.Data.SlopeFilter = SlopeFilter;
output.Data.I_Slope_30 = I_Slope_30;
output.Data.I_Slope_70 = I_Slope_70;
output.Data.dydx = dydx;
output.Data.Median_Trace_Early = zeroed_median_early;
output.Data.Median_Trace_Late = zeroed_median_late;

%% Save Data
Penn_SlopeVal = Penn_SlopeVal';
Max_Amp_Val = Max_Amp_Val';
cd ..
save('All_Data','output')
save('SlopeVals_mV_s','Penn_SlopeVal')
save('AmplitudeVals_mV','Max_Amp_Val')
savefig('fEPSP Graphs')
end
