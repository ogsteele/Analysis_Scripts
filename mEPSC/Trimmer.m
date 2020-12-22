%% Trimmer
% Author: O.G.Steele
% Data: 22.12.20

% Simple script to use ephysIO to load a file, trim it if necessary, then
% use ephysIO to save the file again

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
figure; plot(S.array(:,2))
% trim
[x,~] = ginput(1);
x = floor(x);
a = split(filename,'.');
new_filename = append(char(a(1)),'_trimmed.phy');
close
% save new file
ephysIO(new_filename,S.array(1:x,:),S.xunit,S.yunit,S.names,S.notes,'int16')
% plot the new one
S = ephysIO(new_filename);
figure; plot(S.array(:,2)); title('trimmed trace')
clear