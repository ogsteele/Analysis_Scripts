%% Notes
    % Include a catch statement if there are no filepaths in the .txt
    % Wrap in the script to automatically do the others too
    % Plot EVERYTHING on a master plot (amplitude plots too)
    % Consider integrating the summary.txt stuff too

%% Code 

% PURPOSE:
    % to read through the file paths and compile the events to create new
    % ensemble averages of each of MIXED, AMPAR, NMDAR traces. 

% clear the workspace
clear; close all;

% reference the starting location
home = pwd;

% select file
file = uigetfile('*.txt'); % select file
warning('off') % mute warnings
paths = table2cell(readtable(file)); % read in filepaths
paths = char(paths(:,1)); % convert to character arrays

% Analysis loop
for i = 1:size(paths,1)
    cd(fullfile(strtrim(paths(i,:)),'output_2022')) % cd to the first dir
    % extract AMPAR
    A = ephysIO('AMPAR_ensemble.phy');
    AMPAR(:,i) = A.array(:,2)*1e12; % array in pA
    AMPAR_amp(i) = min(A.array(:,2)); % amp in pA
    % extract Mixed
    M = ephysIO('mixed_ensemble.phy');
    MIXED(:,i) = M.array(:,2)*1e12; % array in pA
    MIXED_amp(i) = min(M.array(:,2)); % amp in pA
    % extract NMDAR (singular, so out of for loop)
    N = ephysIO('NMDAR_ensemble.phy');
    NMDAR(:,i) = N.array(:,2)*1e12; % array in pA
    NMDAR_amp(i) = min(N.array(:,2)); % amp in pA
end

% singular, so out of the loop
time = N.array(:,1); % time in s

% generate averages of the above
AMPAR_ave = median(AMPAR,2);
MIXED_ave = median(MIXED,2);
NMDAR_ave = median(NMDAR,2);

AMPAR_amp_ave = mean(AMPAR_amp);
MIXED_amp_ave = mean(MIXED_amp);
NMDAR_amp_ave = mean(NMDAR_amp);

% plot traces
figure; set(gcf, 'Position',  [100, 100, 1200, 400]);
    % plot the mixed traces 
    subplot(1,3,1); plot(time,MIXED,'color',[0.5 0.5 0.5 0.4],'HandleVisibility','off'); 
    hold on; plot(time,MIXED_ave,'color', 'blue', 'linewidth' ,2)
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)'); legend('Mixed mEPSC','linewidth',1)
    xlim([0 0.1]);ylim([-30 10])

    % plot the AMPAR traces 
    subplot(1,3,2); plot(time,AMPAR,'color',[0.5 0.5 0.5 0.4],'HandleVisibility','off'); 
    hold on; plot(time,AMPAR_ave,'color', 'red', 'linewidth' ,2)
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)'); legend('AMPAR mEPSC','linewidth',1)
    xlim([0 0.1]);ylim([-30 10])
    
    % plot the NMDAR traces 
    subplot(1,3,3); plot(time,NMDAR,'color',[0.5 0.5 0.5 0.4],'HandleVisibility','off'); 
    hold on; plot(time,NMDAR_ave,'color', 'yellow', 'linewidth' ,2)
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)'); legend('NMDAR mEPSC','linewidth',1)
    xlim([0 0.1]);ylim([-30 10])
    
    % create group title
    file = strrep(file,'_',' '); % replace the underscores
    sgtitle(file(1:end-4)) % name the group after the filename
    
    % test run of the subtraction of the two averages
    %figure; plot(time,NMDAR,'color',[0.5 0.5 0.5 0.4],'HandleVisibility','off'); 
    %hold on; plot(time,MIXED_ave-AMPAR_ave,'color', 'yellow', 'linewidth' ,2)
    %box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    %ylabel('Amplitude (pA)'); xlabel('Time (s)'); legend('alt_NMDAR mEPSC','linewidth',1)
    %xlim([0 0.1]);ylim([-30 10])

% return to home
cd(home)