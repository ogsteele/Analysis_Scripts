function [splits, nf_splits] = TP_null_MedianF_function(Param,splits,x,poles)

%% exclude data before median save

%% plan
% split data
% plot wave 
% select exclusion zone
% apply exclusion zone to all waves
% median filter the data without the exclusion zone
% add back in the excluded region
% stitch backinto one whole recording
% convert to function
% integrate into RS_Comp_OGS

%% Parameters

% % warning off
% warning('off','MATLAB:colon:nonIntegerIndex')
% % Amplifier / Sampling settings
% Param.amp_scalef = 0.5; % default scale factor for amplifier used in V / nA
% Param.sample_rate = 20000; % in Hz
% Param.amplifier_scale = 2e09; % actually 2e-09, however 2e09 reverses this
% Param.amplifier_gain = 100; % gain set by experimenter (default is 100)
% Param.split_length = 10; % in seconds, frequency of test pulses
% % Test pulse settings
% Param.pulse_amp = -0.002; % in V
% Param.Vh = -70; % holding voltage in mV
% Param.voltage_step = 2; % in mV
% Param.pulse_duration = 10; % in ms
% Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
% % Compensation settings
% Param.Vrev = 0; % reversal potential in mV (close to zero for AMPAR/NMDAR)
% Param.des_Rs = 8; % desired series resistance for recordings to be compensated to
% %% Load in Data
% % Load in w/ ePhysIO
% % Select raw trace to visualise
% title_str = "1. Select raw file of recording";
% if ~ispc; menu(title_str,'OK'); end
% clear('title_str')
% [file,path,~] = uigetfile('*.*','1. Select raw file of recording');
% % Display file selection selection
% if isequal(file,0)
%    disp('User selected Cancel')
%    % If user selects cancel here, script will end here.
%    return
% else
%     disp(['User selected ', fullfile(path, file)])
%     filename = file;
%     % Navigate to directory and load file with ephysIO
%     cd(path)
%     S = ephysIO(file);
% end
% % tidy workspace
% clear('path','ans')
% 
% %% What gain level was the recording performed at
% % ask how much to trim
% prompt = {'Gain value of the raw recording'};
% dlgtitle = 'Input Gain setting';
% dims = [1 50];
% def = {'100'}; % default gain value in my case is 100
% gain = inputdlg(prompt,dlgtitle,dims,def);
% gain = str2double(gain); % convert to number
% 
% % overwrite the gain value set in parameters
% Param.amplifier_gain = gain;
% %% Split the data into ten second waves
% 
% % convert the data in pA
% S.array(:,2) = (S.array(:,2)/(Param.amp_scalef * Param.amplifier_gain)*1000); % in pA
% % calculate the length of the recording
% length_raw = length(S.array(:,2)); % whole recording in data points
% length_seconds = length_raw/Param.sample_rate; % whole recording in seconds
% n_splits = floor(length_seconds/Param.split_length); % number of splits, rounded down 
% % calculate the number of data points per split
% length_raw_split = Param.split_length * Param.sample_rate;
% 
% % create list of start and end points determined by the sample rate, length
% % of split, length of recording
% % round down n_splits to not deal with incomplete sections and preallocate
% start(:,1) = zeros(n_splits,1);
% finish(:,1) = zeros(n_splits,1);
% % define the starting numbers
% start(1,1) = 1;
% finish(1,1) = 1 + length_raw_split-1;
% for i = 2:n_splits
%     start(i,1) = finish(i-1,1) + 1;
%     finish(i,1) = start(i,1) + length_raw_split-1;
% end
% 
% % create a matrix of the split recording and preallocate
% splits = zeros(length_raw_split,n_splits);
% for s = 1:n_splits
%     splits(:,s) = S.array(start(s):finish(s),2);
% end
% %% select region to exclude
% 
% % get the start time of the test pulse
% figure; plot(splits(:,1))
% title("2. Zoom into the test pulse BEFORE clicking enter to select a single start point")
% pause
% title_str = "2. Zoom into the test pulse BEFORE clicking enter to select a single start point";
% disp("2. Zoom into the test pulse BEFORE clicking enter to select a single start point")
% if ~ispc; menu(title_str,'OK'); end
% clear('title_str')
% [x,y] = ginput(1);
% close

%% exlcudedata, filter and restitch
% % warning off
warning('off','MATLAB:colon:nonIntegerIndex')
nf_splits = splits;

for i = 1:size(splits,2)
    array = splits(:,i);
    time = 0:5.0000e-05:Param.split_length;
    time = time(1:Param.sample_rate*Param.split_length)';

    % select the range to exclude
    % tf = excludedata(ydata,xdata,'range',[100, 150]);
    range_dp = [x(1), x(1)+(0.035*Param.sample_rate)];
    range_s = range_dp/Param.sample_rate;
    
    tf = excludedata(array,time,'range',range_s); % not excl.
    tx = ~excludedata(array,time,'range',range_s); % excl.
    
    % filter the region, bar exclusions
    yf = medianf(array(tf), time(tf) ,poles); % filtered array
    
    % now add back in the bit you didn't want to be filtered
    full = vertcat(yf(1:x(1)),array(tx),yf((x(1)):end));
    splits(:,i) = full;
end

% bug testing plots
% plot region with exclusion zone
% % plot(time(tf),array(tf))
% % plot the excluded region over the top (another color by default
% % hold on; plot(time(tx),array(tx))
% % title('Region excluded from median filtering')
% % figure; plot(full)
% % 
% % figure; plot(tf); hold on; plot(tx); ylim([-0.5 1.5]); hold on
% % legend('tf','tx')
% % xline(range_dp(1)); hold on; xline(range_dp(2)); hold off
% % legend('tf','tx')
% % title('tf / tx')
