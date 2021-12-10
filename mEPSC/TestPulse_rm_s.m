function [tprm,med,x,x2] = TestPulse_rm_s(data)
% standalone version of the test pulse remove function
% input arguments
    % data should be array output from ephysIO
% remove test pulse and plot without

%% Set parameters
prompt = {'Enter Amplifier Scale Factor:',...
    'Enter Sample Rate (Hz):',...
    'Enter Gain Value of Recording:',...
    'Enter Split Length (s):',...
    'Enter Test Pulse Amplitude (V):',...
    'Enter Holding Potential (mV):',...
    'Enter Test Pulse Voltage Step (mV):',...
    'Enter Test Pulse Duration (ms):',...
    'Enter Reversal Potential (mV):',...
    'Enter Desired Rs Value (MOhms):'};
dlg_title = 'Input parameters';
num_lines = 1;
def = {'0.5','20000','100','10','-0.002','-70','2','10','0','8'};
answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
answer = str2double(answer);

% Amplifier / Sampling settings
Param.amp_scalef = answer(1); % default scale factor for amplifier used in V / nA
Param.sample_rate = answer(2); % in Hz
Param.amplifier_gain = answer(3); 
Param.split_length = answer(4); % in seconds, frequency of test pulses
% Test pulse settings
Param.pulse_amp = answer(5); % in V
Param.Vh = answer(6); % holding voltage in mV
Param.voltage_step = answer(7); % in mV
Param.pulse_duration = answer(8); % in ms
Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
% Compensation settings
Param.Vrev = answer(9); % reversal potential in mV (close to zero for AMPAR/NMDAR)
Param.des_Rs = answer(10); % desired series resistance for recordings to be compensated to


%% Split the data into ten second waves

% convert the data in pA
S.array = data;
S.array(:,2) = (S.array(:,2)/(Param.amp_scalef * Param.amplifier_gain)*1000); % in pA
% calculate the length of the recording
length_raw = length(S.array(:,2)); % whole recording in data points
length_seconds = length_raw/Param.sample_rate; % whole recording in seconds
n_splits = floor(length_seconds/Param.split_length); % number of splits, rounded down 
% calculate the number of data points per split
length_raw_split = Param.split_length * Param.sample_rate;

% create list of start and end points determined by the sample rate, length
% of split, length of recording
% round down n_splits to not deal with incomplete sections and preallocate
start(:,1) = zeros(n_splits,1);
finish(:,1) = zeros(n_splits,1);
% define the starting numbers
start(1,1) = 1;
finish(1,1) = 1 + length_raw_split-1;
for i = 2:n_splits
    start(i,1) = finish(i-1,1) + 1;
    finish(i,1) = start(i,1) + length_raw_split-1;
end

% create a matrix of the split recording and preallocate
splits = zeros(length_raw_split,n_splits);
for s = 1:n_splits
    splits(:,s) = S.array(start(s):finish(s),2);
end

% tidy up workspace
clear('start','finish','i','s','length_seconds','length_raw','length_raw_split','n_splits')

% get the start time of the test pulse
figure; plot(splits(:,1))
title("2. Zoom into the test pulse BEFORE clicking enter to select a single start point")
pause
title_str = "2. Zoom into the test pulse BEFORE clicking enter to select a single start point";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[x,~] = ginput(1);
close

% create pulse paramaters from ginput selection
Param.pulse_start = round(x);
Param.pulse_end = Param.pulse_start + Param.pulse_points;                                     % note, this will only cover the downward transiet, not the upward transient
Param.pulse_window = Param.pulse_start:Param.pulse_end;

% calculate the pre pulse mean to fill the gaps
tprm_splits = splits;
base = zeros(size(tprm_splits,2),1);
median_YF = zeros(size(tprm_splits,2),1);
for i = 1:size(tprm_splits,2)
    base(i) = trimmean(tprm_splits(1:Param.pulse_start,i),33,'floor'); % baseline average (in pA)
    tprm_splits(Param.pulse_start:(Param.pulse_end+250),i) = base(i); % swap pulse for base
    median_YF(i) = median(tprm_splits(:,i)); % calculate the median of the splits
end
% concatenate tprm_splits
tprm_splits_conc = vertcat(tprm_splits(:));
% apply a filter to clean the look of the data
t = (0:size(tprm_splits_conc,1)-1)';
t = t./Param.sample_rate;
YF = filter1(tprm_splits_conc, t, 0, 300);
x = linspace(0,(size(tprm_splits_conc,1)/20000),size(tprm_splits_conc,1));
%figure
%plot(x,YF,'Color',[0 0.4470 0.7410 0.2]) % to make lighter
%title('filter')
x2 = linspace(0,(size(tprm_splits_conc,1)/20000),size(median_YF,1));
%hold on
%plot(x2,smooth(median_YF),'linewidth',4,'color',[0 0.4470 0.7410]) % not made lighter
tprm = YF;
med = median_YF;
end