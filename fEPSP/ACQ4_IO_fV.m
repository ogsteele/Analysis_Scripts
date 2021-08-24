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
%path =  'D:\Downloads\OneDrive_2021-08-23\LTP - Base_001'; % whilst in development
path = uigetdir;
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
xlabel('Data Points'); ylabel('Amplitude (V)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
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
hold on; xline(x(1),'color','red','linestyle','--'); hold off;
hold on; xline(x(2),'color','red','linestyle','--'); hold off;
title('fEPSP traces overlaid with points of interest');

% get the minimum value from this region (FV Peak Amp)
for i = 2:size(array,2)
    [max_fv_peak(i-1), max_fv_ind(i-1)] = min(array(x(1):x(2),i));
    max_fv_ind(i-1) = max_fv_ind(i-1) + (x(1) - 1);
    mean_fv_peak(i-1) = mean(array(max_fv_ind(i-1)-1:max_fv_ind(i-1)+1,i));
    hold on; plot(max_fv_ind(i-1),max_fv_peak(i-1),'ro','MarkerSize',10,'linewidth',3); hold off;
end

% FIELD POTENTIAL
% Seperate out window of interest
% make the start of the fEPSP the maximum point AFTER the FV peak
for i = 2:size(array,2)
    [fEPSP_start_val(i-1), fEPSP_start_ind(i-1)] = max(array(max_fv_ind(i-1):x(2),i));
    fEPSP_start_ind(i-1) = fEPSP_start_ind(i-1) + (max_fv_ind(i-1) - 1);
    hold on; plot(fEPSP_start_ind(i-1),fEPSP_start_val(i-1),'go','MarkerSize',10,'linewidth',3); hold off;
end

% calculate the averaged trace from all waves and plot the end
av_trace = (mean(array(:,2:end),2));
fEPSP_end = min(find((abs(av_trace(x(2):end)-0) < 0.00000001))) +x(2);
win_size = fEPSP_end - x(2);
hold on; xline(fEPSP_end,'color','red','linestyle','--'); hold off;

% invert the trace and use findpeaks with constraints
for i = 2:size(array,2)
    fEPSP_Window(:,i-1) = array(fEPSP_start_ind(i-1):(fEPSP_start_ind(i-1) + win_size),i)*-1; % flip the window
    % find the first and largest peak seperated from eachother by 500 and
    % at least 10 points wide to remove noise/spikes
    [pv(i-1),pi_t(i-1)] = findpeaks(fEPSP_Window(:,i-1),'MinPeakDistance',500,'MinPeakWidth',10,'SortStr','descend','NPeaks',1);
    pi(i-1) = pi_t(i-1) + (fEPSP_start_ind(i-1)-1);
    hold on; plot(pi(i-1),pv(i-1)*-1,'yo','MarkerSize',10,'linewidth',3); hold off; 
end

% plot a new window for the slope value
win_time = Time(1:win_size+1); % calculate the time vector for the size of the fEPSP window
fEPSP_Window_f150 = filter1(fEPSP_Window(:,:),win_time,0,150); % 150 Hz Filter
xdiff = fEPSP_Window_f150(1,:);
figure;
for i = 1:size(fEPSP_Window,2)
    xzeroed_fEPSP_Window(:,i) = (fEPSP_Window_f150(:,i) - xdiff(i))*-1;
    hold on; plot(xzeroed_fEPSP_Window(1:pi_t(i),i)); hold off
end
xlabel('Data Points'); ylabel('Amplitude (V)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('fEPSP start to peak amplitude');
xlim([0 300])

% set slope range values
slope_s = 0.20;
slope_e = 0.60;

% calculate the amplitude and index at above values % maximum amplitude along
% the slope and plot as blue and green circles
for i = 1:size(fEPSP_Window,2)
    Amp_s(i) = min(xzeroed_fEPSP_Window(:,i))*slope_s;
    Amp_e(i) = min(xzeroed_fEPSP_Window(:,i))*slope_e;
    Amp_s_ind(i) = min(find((abs(xzeroed_fEPSP_Window(1:pi_t(i),i)-Amp_s(i)) < 0.000005)));
    Amp_e_ind(i) = (min(find((abs(xzeroed_fEPSP_Window(1:pi_t(i),i)-Amp_e(i)) < 0.000005))))+2;
    hold on; plot(Amp_s_ind(i),Amp_s(i),'bo','MarkerSize',10,'linewidth',3); hold off;
    hold on; plot(Amp_e_ind(i),Amp_e(i),'bo','MarkerSize',10,'linewidth',3); hold off;
    hold on; plot((Amp_s_ind(i):1:Amp_e_ind(i)),...
        xzeroed_fEPSP_Window(Amp_s_ind(i):Amp_e_ind(i),i),...
        'linewidth',2,'color','blue'); hold off

    [dydx, ~, ~] = ndiff(xzeroed_fEPSP_Window(Amp_s_ind(i):Amp_e_ind(i),i), win_time(Amp_s_ind(i):Amp_e_ind(i)));
    [Penn_SlopeVal(i), SlopeVal_Index(i)] = min(dydx);
end


% Half Maximal Stimulation Intensity
Stim_Vals_mA = 0:1:9;
pv = pv-(pv(1));
InterpolatedLine = interp1(Stim_Vals_mA,pv,1:0.05:10);
InterpolatedLine = InterpolatedLine - InterpolatedLine(1);
MaxAmp = max(InterpolatedLine);
Half_MaxAmp = MaxAmp/2;
HM_stim_ind = min(find((abs(InterpolatedLine(:)-Half_MaxAmp) < 0.00001)));
Stim_Vals_mA = 0:0.05:9;
HM_Stim = Stim_Vals_mA(HM_stim_ind);
figure;
plot(Stim_Vals_mA,(InterpolatedLine-InterpolatedLine(1)))
xlabel('Stimulation Intensity (mA)'); ylabel('fEPSP Amplitude (V)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('Half Maximal Plot')
yline(Half_MaxAmp)
xline(HM_Stim)

% Annotate with Half Maximal Stimulation Intensity
text = 'Half maximal stimulation intensity (mA) = %f';
str = sprintf(text,HM_Stim);
dim = [0.141071428571428 0.832142857142859 0.43125 0.053571428571429];
annotation('textbox',dim,'String',str,'FitBoxToText','on','linewidth',1,...
    'fontsize',8);

% Stimulation Intensity vs Slope Value
Stim_Vals_mA = 0:1:9;
figure;
plot(Stim_Vals_mA,Penn_SlopeVal*-1);
xlabel('Stimulation Intensity (mA)'); ylabel('Slope Value (V/S)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('Stimulation Intensity vs Slope Value Plot')

% Stimulation Intensity vs Fibre Volley Amplitude
figure;
plot(Stim_Vals_mA,mean_fv_peak*-1);
xlabel('Stimulation Intensity (mA)'); ylabel('Fibre Volley Amplitude (V)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('Stimulation Intensity vs Fibre Volley Amplitude Plot')

% Fibre Volley Amplitude vs Slope Value
figure;
plot(mean_fv_peak*-1,Penn_SlopeVal*-1);
xlabel('Fibre Volley Amplitude (V)'); ylabel('Slope Value (V/S)'); box off; set(gcf,'color','white'); set(gca,'LineWidth',2)
title('Fibre Volley Amplitude vs Slope Value Plot')

%% Define output
cd ..
out.meta.S = S;
out.meta.Meta = Meta;
out.path = path;
out.Stim_Vals_mA = Stim_Vals_mA;
out.FV_Amp = mean_fv_peak*-1;
out.fEPSP_Amp = pv;
out.fEPSP_Slope = Penn_SlopeVal*-1;
out.HM_Stim = HM_Stim;
save('out')
end