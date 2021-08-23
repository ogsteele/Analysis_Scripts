function ACQ4_IO_fV
% Author: O.G. Steele
% Date of Creation: 21.02.19
% Updated: 19.08.21

% New iteration to account for the fibre volley and then plot fibre volley
% to Penn slope ratio input output curves.

%% Params

stim_point = 1000; % number of points at which the stimulus artifact starts
%% Obtain Traces
% Clear windows
close all

% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox\images\imdemos');

% Ask user to confirm or change.
path =  'D:\Downloads\OneDrive_2021-08-23\LTP - Base_001'; % whilst in development
%path = uigetdir;
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

% Collect params/meta via ephysIO
S = ephysIO({[path,'/000/Clamp2.ma'],1}); % loads channel 1 by default
Meta = readMeta([path,'/000/Clamp2.ma']); % reads metadata on the first clamp2.ma file, should be consistent across all
Time = S.array(:,1);
Time_Diff = S.xdiff;

% Change directory to first sweep 
% Note, not first folder as that's the master folder
% PREALLOCATE HERE
for i = 1:numberOfFolders
    cd(char(listOfFolderNames(i,:)));
    Data = h5read('Clamp2.ma','/data');
    Trace(:,i) = Data(:,1);
end

%% Adjust Traces 

% adjust sweeps to zero
basemean = mean(Trace(1:stim_point,:));
adjusted = Trace(:,:) - basemean;

% append time and optionally save as ephysIO format

% apply a 500 Hz LPF
array = [Time (filter1(adjusted(:,:), Time, 0, 500))];
%% Identify Features
% NOTE: Programme this section in line with peakscale.m

% FIBRE VOLLEY
figure; 
plot(array(:,2:end))
xlabel('Data Points'); ylabel('Amplitude (mV)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('1. Locate afferent fibre volley and hit enter')

% allow user to zoom and then hit enter when at region to exit zoom
zoom on
waitfor(gcf, 'CurrentCharacter', char(13))
zoom reset
zoom off

% get points pre and post fibre volley and plot regions
title('2. Isolate the fibre volley with two clicks')
[x,~] = ginput(2);
x = round(x);
hold on; xline(x(1),'color','red','linestyle','--')
hold on; xline(x(2),'color','red','linestyle','--')
% get the minimum value from this region (FV Peak Amp)

for i = 2:size(array,2)
    [max_fv_peak(i-1), max_fv_ind(i-1)] = min(array(x(1):x(2),i-1));
    max_fv_ind(i-1) = max_fv_ind(i-1) + (x(1) - 1);
    mean_fv_peak(i-1) = mean(array(max_fv_ind(i-1)-1:max_fv_ind(i-1)+1,i-1));
    hold on; plot(max_fv_ind(i-1),max_fv_peak(i-1),'ro','MarkerSize',10,'linewidth',3)
end

% FIELD POTENTIAL
% Seperate out window of interest
% make the start of the fEPSP the maximum point AFTER the FV peak
for i = 2:size(array,2)
    [fEPSP_start_val(i-1), fEPSP_start_ind(i-1)] = max(array(max_fv_ind(i-1):x(2),i-1));
    fEPSP_start_ind(i-1) = fEPSP_start_ind(i-1) + (max_fv_ind(i-1) - 1);
    hold on; plot(fEPSP_start_ind(i-1),fEPSP_start_val(i-1),'go','MarkerSize',10,'linewidth',3)
end

% calculate the averaged trace from all waves and plot the end
av_trace = (mean(array(:,2:end),2));
fEPSP_end = min(find((abs(av_trace(x(2):end)-0) < 0.00000001))) +x(2);
win_size = fEPSP_end - x(2);
hold on; plot(fEPSP_end,0,'bo','MarkerSize',10,'linewidth',3);

% find the maximum amplitude of each trace
for i = 2:size(array,2)
    extra_filtered_array(:,i-1) = array(fEPSP_start_ind(i-1):(fEPSP_start_ind(i-1) + win_size),i);
end
extra_filtered_array_time = Time(1:win_size+1);

for i = 2:size(array,2)
    [fEPSP_max_amp_val(i-1),fEPSP_max_amp_ind(i-1)] = min(array(fEPSP_start_ind(i-1):fEPSP_end,i));
    fEPSP_max_amp_ind(i-1) = fEPSP_max_amp_ind(i-1) + (fEPSP_start_ind(i-1)-1);
    hold on; plot(fEPSP_max_amp_ind(i-1),fEPSP_max_amp_val(i-1),'yo','MarkerSize',10,'linewidth',3); 
end


fEPSP_Window = array((fEPSP_start:fEPSP_end),2:end);
num_points = size(fEPSP_Window,1);
fEPSP_Window_Time = linspace(0,Time_Diff * num_points,num_points);

% assume that the field potential begins immediately after the fibre volley

%% Get Points
% plot the sweeps and get points pre and post field potential
plot(adjusted(stim_point:5500,:));
ylim([-0.002 0.0005]);
xlim([0 2000]);
hold on
[x,y] = getpts;
close

% Seperate out window of interest
fEPSP_start = round(x(1)+stim_point);
fEPSP_end = round(x(2)+stim_point);
fEPSP_Window = adjusted((fEPSP_start:fEPSP_end),:);

% Calculate time of window
%Time_Conv = length(fEPSP_Window)/samplerate;
%xtime = linspace(0,Time_Conv,length(fEPSP_Window));
xtime = 0:Time_Diff:Time_Diff*(length(fEPSP_Window)-1);

% Plot fEPSP Traces
%subplot(1,2,1)
figure
for l = 1:10
    plot(xtime,(smooth(smooth((fEPSP_Window(:,l)*1000)))))
    hold on
end

ylabel('Amplitude (mV)','fontsize',16);
xlabel('Time (ms)','fontsize',16);
%xlim([0 Time_Conv]);
title('fEPSP I/V Traces','fontsize',18)
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)
legend({'0 V','10 V','20 V','30 V','40 V','50 V','60 V','70 V','80 V',...
    '90 V'},'Location','southeast','linewidth',1,'fontsize',14);
hold off


% Index of maximum amplitude, and incidentally the peak amplitude of each
% stimulation intensity
for o = 1:10
    [M(o),I(o)] = min((fEPSP_Window(:,o)));
end
M = M*1000;

% index of 30% and 70% of maximum amplitudes
I_slope_30 = round(I*0.3);
I_slope_70 = round(I*0.7);

% calculate fEPSP slope between 30% and 70% peak amplitude
for k = 1:10
    slope(k) = mean(gradient(fEPSP_Window((I_slope_30(k):I_slope_70(k)),k)));
end

% Half Maximal Stimulation Intensity
StimIntensity = linspace(0,90,10)';
HM = interp1(M,StimIntensity,(min(M)/2));
MaxAmp = min(M);

% Plot I/O Curve
% subplot(1,2,2)
figure
plot(StimIntensity,M');
hold on
ylabel('Amplitude (mV)','fontsize',16);
xlabel('Stimulation Intensity (V)','fontsize',16);
title('I/O Curve','fontsize',18);
box off
set(gcf,'color','w')
set(gca,'linewidth',2,'fontsize',14)

% Annotate with Half Maximal Stimulation Intensity
text = 'Half maximal stimulation intensity (V) = %f';
str = sprintf(text,HM);
dim = [0.441071428571428 0.832142857142859 0.43125 0.053571428571429];
annotation('textbox',dim,'String',str,'FitBoxToText','on','linewidth',1,...
    'fontsize',14);

% Create output structure
output.MaxAmp = MaxAmp;
output.HalfMax = HM;
end