function [LTP_Data] = ACQ4_LTP

% Author: O.G. Steele
% Date of Creation: 21.02.19
% Updated: 28.02.19

%% Syntax

%% Update Notes

% 28.02.19
% Data contained in structure, saved at site level as only one LTP
% experiment should be carried out per slice. 


%% Define Sampling rate in kHz
samplerate = 40;

%% Clear workspace, close figures
clear all
close all

%% PRE STIMULUS
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
adjusted1 = Trace(:,:) - basemean;

plot(adjusted1(3500:5500,:));
ylim([-0.004 0.005]);
xlim([0 2000]);
hold on
[x,y] = getpts;
close

% Seperate out window of interest
fEPSP_start = round(x(1)+3500);
fEPSP_end = round(x(2)+3500);
fEPSP_Window1 = adjusted1((fEPSP_start:fEPSP_end),:);
[H_fw,numofrecordings] = size(fEPSP_Window1);

% Index of maximum amplitude, and incidentally the peak amplitude of each
% stimulation intensity
for o = 1:numofrecordings
    [M1(o),I(o)] = min((fEPSP_Window1(:,o)));
end
M1 = M1*1000;

% index of 30% and 70% of maximum amplitudes
I_slope_30 = round(I*0.3);
I_slope_70 = round(I*0.7);

% calculate fEPSP slope between 30% and 70% peak amplitude
for k = 1:numofrecordings
    pre_slope(k) = mean(gradient(smooth(smooth(smooth(fEPSP_Window1((I_slope_30(k):I_slope_70(k)),k))))));
end

%% POST STIMULUS

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
adjusted2 = Trace(:,:) - basemean;

plot(adjusted2(3500:5500,:));
ylim([-0.004 0.005]);
xlim([0 2000]);
hold on
[x,y] = getpts;
close

% Seperate out window of interest
fEPSP_start = round(x(1)+3500);
fEPSP_end = round(x(2)+3500);
fEPSP_Window2 = adjusted2((fEPSP_start:fEPSP_end),:);
[H_fw,numofrecordings] = size(fEPSP_Window2);

% Index of maximum amplitude, and incidentally the peak amplitude of each
% stimulation intensity
for o = 1:numofrecordings
    [M2(o),I(o)] = min((fEPSP_Window2(:,o)));
end
M2 = M2*1000;

% index of 30% and 70% of maximum amplitudes
I_slope_30 = round(I*0.3);
I_slope_70 = round(I*0.7);

% calculate fEPSP slope between 30% and 70% peak amplitude
for k = 1:numofrecordings
    post_slope(k) = mean(gradient(smooth(smooth(smooth(fEPSP_Window2((I_slope_30(k):I_slope_70(k)),k))))));
end

%% Plotting LTP

% concatenate slop values pre and post stimulus
slope_vals = [pre_slope post_slope];
time = linspace(-10,60,length(slope_vals));

% normalise
slope_basemean = mean(slope_vals(1:10));
norm_percent_slope_vals = (slope_vals/slope_basemean)*100

% plot slope vals
plot(time,norm_percent_slope_vals);
xlabel('Time (mins)');
ylabel('fEPSP gradient % of baseline')

%% % Change directory up two
mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end-1)); 

%% Assign Variables to Structure
LTP_Data.Date = mydir(idcs(9):idcs(10));
LTP_Data.Slice = mydir(idcs(10):idcs(11));
LTP_Data.Site = mydir(idcs(11):idcs(12));
LTP_Data.AdjustedRawData_preTS = adjusted1;
LTP_Data.AdjustedRawData_postTS = adjusted2;
LTP_Data.fEPSP_Window_preTS = fEPSP_Window1;
LTP_Data.fEPSP_Window_postTS = fEPSP_Window2;
LTP_Data.PeakAmp_preTS = M1;
LTP_Data.PeakAmp_postTS = M2;
LTP_Data.SlopeVals_preTS = pre_slope;
LTP_Data.SlopeVals_postTS = post_slope;
LTP_Data.SlopeVals_all = slope_vals;

%% Saving Structure

% Change directory up two 
cd(newdir)
save LTP_Data
clearvars -except LTP_Data

