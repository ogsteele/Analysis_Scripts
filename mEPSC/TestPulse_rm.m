%% TestPulse_rm

%% Parameters
% Amplifier / Sampling settings
Param.sample_rate = 20000; % in Hz
Param.amplifier_scale = 2e09; % actually 2e-09, however 2e09 reverses this
Param.amplifier_gain = 100;
Param.split_length = 10; % in seconds, frequency of test pulses
% Test pulse settings
Param.pulse_amp = -0.002; % in V
Param.Vh = -45; % holding voltage in mV
Param.commander_scale = 20; % in mV/V
Param.voltage_step = 2; % in mV
Param.pulse_duration = 10; % in ms
Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
Param.pulse_start = 19100;
Param.pulse_end = Param.pulse_start + Param.pulse_points;                                     % note, this will only cover the downward transiet, not the upward transient
Param.pulse_window = Param.pulse_start:Param.pulse_end;
% Compensation settings
Param.Vrev = 0; % reversal potential in mV (close to zero for AMPAR/NMDAR)
Param.des_Rs = 8; % desired series resistance for recordings to be compensated to

%% Data Selection
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


% For use during development - replacing the file dialog
%cd('/Users/ogsteele/Library/Mobile Documents/com~apple~CloudDocs/DPhil/Analysis/Analysis_Scripts/mEPSC')
%S = ephysIO('20201021_000.tdms');
%% Split the data into ten second waves
% convert the data in pA (like seen in nidaq_scope)
S.array(:,2) = S.array(:,2) * (Param.amplifier_scale * Param.amplifier_gain);
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

%% remove test pulse and plot without
% calculate the pre pulse mean to fill the gaps
tprm_splits = splits;
base = zeros(size(splits,2),1);
for i = 1:size(splits,2)
    base(i) = trimmean(tprm_splits(1:Param.pulse_start,i),33,'floor'); % baseline average (in pA)
    tprm_splits(Param.pulse_start:(Param.pulse_end+250),i) = base(i); % swap pulse for base
    median_YF(i) = median(tprm_splits(:,i)); % calculate the median of the splits
end
% concatenate tprm_splits
tprm_splits = vertcat(tprm_splits(:));
figure
plot(tprm_splits)
title('raw')
% apply a filter to clean the look of the data
t= (0:size(tprm_splits,1)-1)';
t = t./Param.sample_rate;
YF = filter1(tprm_splits, t, 0, 1000);
figure
plot(YF,'Color',[0 0.4470 0.7410 0.2]) % make lighter
title('filter')
% interpolate median trace
x = linspace(0,size(tprm_splits,1),180);
hold on
plot(x,smooth(median_YF),'linewidth',4,'color',[0 0.4470 0.7410])
