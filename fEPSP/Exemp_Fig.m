%% Exemplary fEPSP Plots
% Simple script to plot median fEPSP traces for figure

%% Parameters
% Define Sample Rate in kHz
samplerate = 40;


%% Analysis
% Clear windows
close all

% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox\images\imdemos');

% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end

% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);

% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ':');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

% Change directory to first sweep 
% Note, not first folder as that's the master folder
for i = 2:numberOfFolders
    cd(char(listOfFolderNames(:,i)));
    Data = h5read('Clamp2.ma','/data');
    Trace(:,i-1) = Data(:,1);
end

% adjust sweeps to zero and smooth to compensate for noise
basemean = mean(Trace(1:3500,:));
adjusted = Trace(:,:) - basemean;

% plot the sweeps and get points pre and post field potential
plot(adjusted(3500:5500,:));
ylim([-0.004 0.005]);
xlim([0 2000]);
hold on
[x,y] = getpts;
close

% Seperate out window of interest
fEPSP_start = round(x(1)+3500);
fEPSP_end = round(x(2)+3500);
fEPSP_Window = adjusted((fEPSP_start:fEPSP_end),:);
fEPSP_Window = fEPSP_Window * 1e3; % convert to mV

% Calculate time of window
Time_Conv = length(fEPSP_Window)/samplerate;
xtime = linspace(0,Time_Conv,length(fEPSP_Window));

% Plot fEPSP Traces
figure
for l = 1:size(fEPSP_Window,2)
    plot(xtime,(smooth(smooth((fEPSP_Window(:,l)*1000)))))
    hold on
end

ylabel('Amplitude (mV)','fontsize',16);
xlabel('Time (ms)','fontsize',16);
xlim([0 Time_Conv]);
title('fEPSP I/V Traces','fontsize',18)
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
hold off

%% Shaded error bar and median trace
% Plot Median Trace
% Notes: Plots median of the first 6th, middle 6th and final 6th
% In the case of a 60 minute recording, 1/6th = 10 minutes = 20 sweeps
figure
Time_Diff = 1/(samplerate*1e3);
num_points = size(fEPSP_Window,1);
fEPSP_Window_Time = linspace(0,Time_Diff * num_points,num_points);
% define the early traces (first 20 sweeps)
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
ylabel('Amplitude (mV)')
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)