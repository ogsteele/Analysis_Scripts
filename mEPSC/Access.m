%% Access Plot
% Author: O.G. Steele
% Date: 09.10.20
% Description: 
%       Simple MATLAB Script to reveal the shift in access resistance

%% Load in Data

% Select baseline.txt file 
title_str = "1. Select access.txt for the recording of choice";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[file,path,~] = uigetfile('*.txt');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    cd(path)
    % Load the baseline file and seperate out the values within
    Baseline = dlmread(file,' ');
    Time_s = Baseline(:,1);
    Access_mO = Baseline(:,2)*(10^-6);
    % Plot the baseline and calculate difference as the mean of the first
    % third, compared to the mean of the last third
    access_plot = figure;
    plot(Time_s,Access_mO,'Linewidth',2)
    early_access = mean(Access_mO(1:(round(length(Access_mO)/3))));
    late_access = mean(Access_mO((round(length(Access_mO)/3)*2:end)));
    Shift = (late_access-early_access);
    dim = [.2 .5 .3 .3];
    % Stylise graph
    box off
    set(gca,'linewidth',2)
    set(gcf,'color','w');
    % Save the read textfile as a .mat
    matfile = strcat('mat_',file(1:end-4));
    % Calculate the percentage change in shift
    shift_percent = (Shift / early_access)*100;
    access.(matfile).Time_s = Time_s;
    access.(matfile).Access_mO = Access_mO;
    access.(matfile).Shift_mO = Shift;
    access.(matfile).Shift_percent = shift_percent;
    access.(matfile).Shift_vals = [early_access,late_access];
    mkdir Access_Shift_Results
    save('Access_Shift_Results/access.mat','-struct','access');
    annotation('textbox',dim,'String', strcat(num2str(access.(matfile).Shift_percent),' % shift'),'FitBoxToText','on')
    saveas(access_plot,'Access_Shift_Results/shift_fig.pdf')
    % Tidy up
    clearvars -except 'access' 'shift_percent'
end