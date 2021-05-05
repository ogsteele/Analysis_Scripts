%% Clipper

% Clipper clips the beginning, trimmer trims the end ... don't ask
% combine at a later date
%% How much do you want to trim?
% ask how much to trim
prompt = {'Amount of data to trim at the start (s)'};
dlgtitle = 'Input';
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims);
answer = str2double(answer); % convert to number

%% Load in w/ ePhysIO
% Select raw trace to visualise
title_str = "1. Select raw file of recording that requires trimming";
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

%% Plot data, trim and save
% plot
figure; plot(S.array(:,2)), title('raw trace')

% trim late 
a = split(filename,'.');
new_filename = append(char(a(1)),'.phy');
% save new file
ephysIO(new_filename,S.array(answer*20000:end,:),S.xunit,S.yunit,S.names,S.notes,'int16')
% plot the new one
S = ephysIO(new_filename);
figure; plot(S.array(:,2)); title('trimmed trace')

clear

