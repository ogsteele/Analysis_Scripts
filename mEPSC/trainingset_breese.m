%% Training Set Creator
% Simple script, bundled with function, to create a training set.
% Note add in saving of file name and chosen waves at the end
% Also change where the question is for sample rate and split length ...

% ask if user is ready
answer = menu('Ready to create a training set?','nope','yep');
% if answer is no ...
if answer == 1
    % ... end loop and tell the user to type when ready
    disp('Type "trainingset_breese" when ready (if on path)')
    return
% but if the answer is yes ...
elseif answer == 2
    % ... let's go!
    disp('Lets go!')
    % run the extractr function below, with the output as waves
    [waves, sample_rate, filename, wavenum] = extractr;
    % infinite for loop for more waves added onto the end
    for k = 1:50000
        % ask the user if they want to add another
        answer = menu('add another?','nope','yep');
        % if answer is no
        if answer == 1
            % run ephysIO and concatenate time and waves to save the data
            waves = waves(1:(size(waves,1)-1),:); % temp fix for the jump scaling due to Vm change
            time = (1/sample_rate)*[1:size(waves,1)]';
            ephysIO('trainingset.phy',[time (waves*0.01)],'s','A') % .tdms scaling for James
            %ephysIO('trainingset.phy',[time waves],'s','A')
            % display that you're done here
            disp('Training set concatenated. Training set saved in current directory (below)')
            disp(pwd)
            return
        % if the answer is yes though, run it again and go again    
        elseif answer == 2
            [some_vector,~,~,~] = extractr;
            % concatenate next to 
            waves = [waves some_vector];
        end
    end
end


function [waves, sample_rate_Hz, filename, wavenum] = extractr
%% Outputs
% waves = waves of interest
% sample_rate_Hz = sampling rate in Hz
% filename = filename of the recording
% wavenum = wave numbers chosen

%% Load in w/ ePhysIO
% Select raw trace to visualise
title_str = "1. Select raw file of recording";
menu(title_str,'OK');
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
    channel = menu('Select Channel (1 is default)','1','2');
    S = ephysIO({file,channel});
end

%% Select waves
% note: wave 5 has seconds 40 through 50 and sampled at 20 kHz
% select the waves of interest
prompt = {'Wave Number 1', 'Wave Number 2'};
dlgtitle = 'Wave Selection';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
wavenum = answer;
% convert the wave numbers to the start and end times in s and dp
split_length_s = 60;
sample_rate_Hz = 25000;
wave_start_s = (str2double(answer)-1)*split_length_s; 
wave_end_s = wave_start_s+split_length_s;
wave_start_dp = wave_start_s * sample_rate_Hz;
wave_end_dp = wave_end_s * sample_rate_Hz;
% extract the waves of interest
wave_1 = S.array(wave_start_dp(1):wave_end_dp(1),2);
wave_2 = S.array(wave_start_dp(2):wave_end_dp(2),2);
waves = [wave_1 wave_2];
end
