function out = earlypull

%% Exemplary fEPSP Plots
% Simple function to pull the early section of fEPSP recordings and zeroed
% by trace

%% Parameters
% Define Sample Rate in kHz
samplerate = 40;


%% Analysis
% Clear windows
%close all

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
out = adjusted(:,1:20);
figure
plot(out)
