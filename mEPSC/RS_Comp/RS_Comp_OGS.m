%% Rs_Comp_Steele
% Script to load in ephys data, and compensate the series resistance
% changes observed during whole cell patch clamp recordings

% Note, change parameters as relevent

% Dependencies - ensure in path
% ephysIO (https://github.com/acp29/eventer/blob/master/base/ephysIO.m)

%% TO DO 
% legacy filetype in filename, split by delimiter
% optional low pass filter inclusion - not priority
% ensure continuity present with if statements

%% Parameters
Param.sample_rate = 20000; % in Hz
Param.amplifier_scale = 2e09;
Param.amplifier_gain = 100;
Param.split_length = 10; % in seconds, frequency of test pulses
Param.voltage_step = 10; % in mV
%Param.lp_filter_preference = false; % true / false
%Param.lp_filter_cutoff = 1000; % in Hz

%% Load in w/ ePhysIO 
% plus UI for script

% Select raw trace to visualise
title_str = "1. Select raw file of recording";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Select raw file of recording');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    filename = file;
    % Navigate to directory and load file with ephysIO
    cd(path)
    S = ephysIO(file);
    % save raw recording following the ephysIO guidance
    % To save a column-major XY array of electrophysiology data:
    %    ephysIO (filename,array,xunit,yunit,names,notes,datatype)
    % ephysIO((append('uncomp_',filename,'.phy')),S.array,S.xunit,S.yunit,S.names,S.notes,'int16')
    
end
% tidy workspace
clear('file', 'path','ans')

%% Split the data into ten second waves
% calculate the length of the recording
length_raw = length(S.array(:,2));
length_seconds = length_raw/Param.sample_rate;
n_splits = length_seconds/Param.split_length;
% number of data points per split
length_raw_split = length_raw/n_splits;

% create list of start and end points determined by the sample rate, length
% of split, length of recording
% round down n_splits to not deal with incomplete sections and preallocate
start(:,1) = zeros(floor(n_splits),1);
finish(:,1) = zeros(floor(n_splits),1);
% define the starting numbers
start(1,1) = 1;
finish(1,1) = 1 + length_raw_split-1;
for i = 2:floor(n_splits)
    start(i,1) = finish(i-1,1) + 1;
    finish(i,1) = start(i,1) + length_raw_split-1;
end
    
% create a matrix of the split recording and preallocate
splits = zeros(length_raw_split,floor(n_splits));
for s = 1:floor(n_splits)
    splits(:,s) = S.array(start(s):finish(s),2);
end

% tidy up workspace
clear('start','finish','i','s','n_splits','length_seconds','length_raw', ...
    'length_raw_split')

%% Generate necessary whole cell paramaters
