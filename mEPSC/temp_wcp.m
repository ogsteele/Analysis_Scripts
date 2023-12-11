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
end
% tidy workspace
clear('path','ans')

%% Open notes.txt file, if present
% if notes.txt associated with recording exists, print in command window
a = split(file,'.');
notes_file = append(char(a(1)),'_notes.txt');
if exist(notes_file,'file')
    type(notes_file)
    notes_import = readcell(notes_file); % read in notes file
    % check if the notes file being loaded is compatible
    if isequal(char(notes_import(1,1)),'Scale factor to generate data from tdms file')
        params_set = 1;
            % Amplifier / Sampling settings
        Param.scaling = cell2mat(notes_import(1,2)); % scaling to generate data
        Param.amplifier = char(notes_import(2,2)); % Amplifier used
        Param.amp_scalef = cell2mat(notes_import(3,2)); % default scale factor for amplifier used in V / nA 
            temp_sr = char(notes_import(6,2));
        Param.sample_rate = str2double(temp_sr(1:end-3)); % in Hz
            clear temp_sr
        Param.amplifier_gain = cell2mat(notes_import(4,2));
            temp_sl = char(notes_import(10,2));
        Param.split_length = str2double(temp_sl(1:end-2)); % in seconds, frequency of test pulses
            clear temp_sl
        Param.amp_lpf = char(notes_import(5,2)); % low pass filter cut off in Hz
        Param.yunit = char(notes_import(7,2));
        Param.xunit = 'S';
        % Test pulse settings
            temp_Vh = char(notes_import(9,2));
        Param.Vh = str2double(temp_Vh(1:end-11)); % holding voltage in mV
            clear temp_Vh
            temp_vs = char(notes_import(13,2));
        Param.voltage_step = abs(str2double(temp_vs(1:end-11))); % in mV
            clear temp_vs;
            temp_pd = char(notes_import(12,2));
        Param.pulse_duration = (str2double(temp_pd(1:end-2))); % in s
            clear temp_pd
        Param.pulse_points = (Param.pulse_duration*Param.sample_rate);               % convert from ms to data points
        % Compensation settings
        Param.Vrev = 0; % reversal potential in mV (close to zero for AMPAR/NMDAR)
        Param.des_Rs = 8; % desired series resistance for recordings to be compensated to
    else % if notes file is not compatible, enter manually
        params_set = 0;
        disp('Notes file not compatible, check notes file or enter manually')
        prompt = {'Enter Scaling to generate data',...
            'Enter Amplifier Scale Factor:',...
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
        def = {'2e-11','2e-09','25000','100','10','-0.002','-60','2','10','0','8'};
        answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
        answer = str2double(answer);
    end
else
    params_set = 0;
    disp('No notes file present')
    prompt = {'Enter Scaling to generate data',...
        'Enter Amplifier Scale Factor:',...
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
    def = {'2e-11','2e-09','25000','100','10','-0.002','-60','2','10','0','8'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
    answer = str2double(answer);
end

% Set parameters
% prompt = {'Enter Scaling to generate data',...
%     'Enter Amplifier Scale Factor:',...
%     'Enter Sample Rate (Hz):',...
%     'Enter Gain Value of Recording:',...
%     'Enter Split Length (s):',...
%     'Enter Test Pulse Amplitude (V):',...
%     'Enter Holding Potential (mV):',...
%     'Enter Test Pulse Voltage Step (mV):',...
%     'Enter Test Pulse Duration (ms):',...
%     'Enter Reversal Potential (mV):',...
%     'Enter Desired Rs Value (MOhms):'};
% dlg_title = 'Input parameters';
% num_lines = 1;
% def = {'2e-11','2e-09','25000','100','10','-0.002','-60','2','10','0','8'};
% answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
% answer = str2double(answer);

% Amplifier / Sampling settings
if params_set == 0 % if params have not been set as above
    Param.scaling = answer(1); % scaling to generate data
    Param.amp_scalef = answer(2); % default scale factor for amplifier used in V / nA
    Param.sample_rate = answer(3); % in Hz
    Param.amplifier_gain = answer(4); 
    Param.split_length = answer(5); % in seconds, frequency of test pulses
    % Test pulse settings
    Param.Vh = answer(7); % holding voltage in mV
    Param.voltage_step = answer(8); % in mV
    Param.pulse_duration = answer(9); % in ms
    Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
    % Compensation settings
    Param.Vrev = answer(10); % reversal potential in mV (close to zero for AMPAR/NMDAR)
    Param.des_Rs = answer(11); % desired series resistance for recordings to be compensated to
else
end

%% Split the data into ten second waves

% convert the data in pA
S.array(:,2) = (S.array(:,2)*Param.scaling)* 1e12; % in pA


% calculate the length of the recording
length_seconds = length(S.array(:,2))/Param.sample_rate; % whole recording in seconds
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
for j = 1:n_splits
    splits(:,j) = S.array(start(j):finish(j),2);
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
[x,y] = ginput(1);
close

% create pulse paramaters from ginput selection
Param.pulse_start = round(x);
Param.pulse_end = Param.pulse_start + Param.pulse_points; % note, this will only cover the downward transiet, not the upward transient
Param.pulse_window = Param.pulse_start:Param.pulse_end;
pulsewin = S.array(Param.pulse_window,2);

% positive or negative pulse
% Was the recording exposed to LPS/PBS?
%dlgTitle    = 'Upward or downward deflection';
%dlgQuestion = 'Which deflection, upward or downward, was selected?';
%if istep_run == 0
%   condition = questdlg(dlgQuestion,dlgTitle,'LPS (Test)','PBS (Control)','PBS (Control');
%elseif istep_run == 1
%   condition = questdlg(dlgQuestion,dlgTitle,'LPS (Test)','PBS (Control)',(output.condition));
%end

if pulsewin(1) > pulsewin(end)
    pulsedir = 'Negative';
elseif abs(pulsewin(1)) > abs(pulsewin(end))
    pulsedir = 'Positive';
end


% plot first and last window used to calculate wcp parameters
% figure; plot(splits(Param.pulse_window,1))
% hold on
% plot(splits(Param.pulse_window,size(splits,2)))
% legend("first pulse","last pulse")
% title("first and last raw test pulses overlaid")

%% Generate necessary whole cell paramaters
wcp_raw = WCP(splits,Param);
disp('Generating raw whole cell properties ... ')

