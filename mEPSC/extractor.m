function waves = extractor
%% Load in w/ ePhysIO
% Select raw trace to visualise
title_str = "1. Select raw compensated recording";
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
    S = ephysIO(file);
end
% tidy workspace
clear('path','ans')

%% Select waves
% note: wave 5 has seconds 40 through 50 and sampled at 20 kHz
% select the waves of interest
prompt = {'Wave Number 1', 'Wave Number 2'};
dlgtitle = 'Wave Selection';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
% convert the wave numbers to the start and end times in s and dp
wave_start_s = (str2double(answer)-1)*10; 
wave_end_s = wave_start_s+10;
wave_start_dp = wave_start_s * 20000;
wave_end_dp = wave_end_s * 20000;
% extract the waves of interest
wave_1 = S.array(wave_start_dp(1):wave_end_dp(1),2);
wave_2 = S.array(wave_start_dp(2):wave_end_dp(2),2);
waves = [wave_1 wave_2];

% your_result = [];
% for ii = whatever
%  some_vector = some_function(of_something);
% concatenate next to 
%  your_result = [your_result {some_vector}];
% end
% title_str = "Add another recording?";
% one_more = menu(title_str,'Yes','No');
% 1 = Yes, 2 = No
% clear('title_str')
end
