%% Notes
    % Include a catch statement if there are no filepaths in the .txt
    % Wrap in the script to automatically do the others too
    % Plot EVERYTHING on a master plot (amplitude plots too)

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
    % pull out event ensemble information
    cd(fullfile(strtrim(paths(i,:)),'output_2022')) % cd to the first dir
    % extract AMPAR
    A = ephysIO('AMPAR_ensemble.phy');
    AMPAR(:,i) = A.array(:,2)*1e12; % array in pA
    % extract Mixed
    M = ephysIO('mixed_ensemble.phy');
    MIXED(:,i) = M.array(:,2)*1e12; % array in pA
    % extract NMDAR (singular, so out of for loop)
    N = ephysIO('NMDAR_ensemble.phy');
    NMDAR(:,i) = N.array(:,2)*1e12; % array in pA
    NMDAR_amp(i,1) = min(movmean(NMDAR(:,i),100)); % calculate NMDAR amplitude here (in pA)
    %NMDAR amplitude taken as a 100 point moving average of the trace to
    %account for any noise (could filter, but this is quicker)
    
    % pull out event frequency information
    cd ..;
    before_summary = table2cell(readtable(...
        'mlm_before\eventer.output\ALL_events\summary.txt'));
    MIXED_amp(i,1) = cell2mat(before_summary(6,2)); % amplitude in pA
    MIXED_Hz(i,1) = cell2mat(before_summary(4,2)); % frequency in Hz
    MIXED_rise(i,1) = cell2mat(before_summary(8,2)); % rise time constant in ms
    MIXED_decay(i,1) = cell2mat(before_summary(9,2)); % decay time constant in ms
    
    after_summary = table2cell(readtable(...
        'mlm_after\eventer.output\ALL_events\summary.txt'));
    AMPAR_amp(i,1) = cell2mat(after_summary(6,2)); % amplitude in pA
    AMPAR_Hz(i,1) = cell2mat(after_summary(4,2)); % frequency in Hz
    AMPAR_rise(i,1) = cell2mat(after_summary(8,2)); % rise time constant in ms
    AMPAR_decay(i,1) = cell2mat(after_summary(9,2)); % decay time constant in ms
    
end

% singular, so out of the loop
time = N.array(:,1); % time in s

% generate averages of the above traces
AMPAR_ave = median(AMPAR,2);
MIXED_ave = median(MIXED,2);
NMDAR_ave = median(NMDAR,2);

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

% plot the summary data (boxplots)
figure;
    % plot the amplitude data
    subplot(2,2,1)
    Y = [MIXED_amp AMPAR_amp abs(NMDAR_amp)];
    b = boxplot(Y,'labels',{'MIXED','AMPAR','NMDAR'}); box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','blue');
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','red');
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','yellow');
    xlabel('Ensemble composition');ylabel('Amplitude (pA)');
    
    % plot the frequency data
    subplot(2,2,2)
    Y = [MIXED_Hz AMPAR_Hz];
    b = boxplot(Y,'labels',{'MIXED','AMPAR'}); box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:2,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','blue');
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','red');
    xlabel('Ensemble composition');ylabel('Frequency (Hz)');
    
    % plot the rise time data
    subplot(2,2,3)
    Y = [MIXED_rise AMPAR_rise];
    b = boxplot(Y,'labels',{'MIXED','AMPAR'}); box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:2,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','blue');
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','red');
    xlabel('Ensemble composition');ylabel('Rise Time (ms)');
    
    % plot the decay time constant data
    subplot(2,2,4)
    Y = [MIXED_decay AMPAR_decay];
    b = boxplot(Y,'labels',{'MIXED','AMPAR'}); box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:2,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','blue');
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','red');
    xlabel('Ensemble composition');ylabel('Decay time constant (ms)');
    
    % create group title
    file = strrep(file,'_',' '); % replace the underscores
    sgtitle(file(1:end-4)) % name the group after the filename
% return to home
cd(home)