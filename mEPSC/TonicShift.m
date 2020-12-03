%% Tonic Current Shift
% Author: O.G. Steele
% Date: 09.10.20
% Description: 
%       Simple MATLAB Script to reveal the holding current shift after
%       applying NMDAR antagonist during subtraction experiments

%% Load in Data

% Select baseline.txt file 
title_str = "1. Select Baseline.txt for the recording of choice";
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
    Holding_pA = Baseline(:,2)*(10^12);
    % Plot the baseline and calculate difference as the mean of the first
    % third, compared to the mean of the last third
    shift_plot = figure;
    plot(Time_s,Holding_pA,'Linewidth',2)
    early_current = mean(Holding_pA(1:(length(Holding_pA)/3)));
    late_current = mean(Holding_pA((length(Holding_pA)/3)*2:end));
    Shift = (early_current - late_current);
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String', strcat(num2str(Shift),' pA shift'),'FitBoxToText','on')
    % Stylise graph
    box off
    set(gca,'linewidth',2)
    set(gcf,'color','w');
    xlabel('Time (seconds)')
    ylabel('Holding Current (pA)')
    % Save the read textfile as a .mat
    matfile = strcat('mat_',file(1:end-4));
    baseline.(matfile).Time_s = Time_s;
    baseline.(matfile).Holding_pA = Holding_pA;
    baseline.(matfile).Shift_pA = Shift;
    baseline.(matfile).Shift_percent = (Shift*-1 / early_current*-1)*100;
    baseline.(matfile).Shift_vals = [early_current,late_current];
    mkdir Tonic_Shift_Results
    save('Tonic_Shift_Results/baseline.mat','-struct','baseline');
    saveas(shift_plot,'Tonic_Shift_Results/shift_fig.pdf')
    % Tidy up
    clearvars -except 'baseline'
end
    

%Note don't forget about T =
%table2array(readtable('20201021_000_access.txt')) to get text into arrays
    
    