%% Notes
    % Include a catch statement if there are no filepaths in the .txt
    % Wrap in the script to automatically do the others too
    % Plot EVERYTHING on a master plot (amplitude plots too)

%% Code 

% PURPOSE:
    % to read through the file paths and compile the events to create new
    % ensemble averages of each of MIXED, AMPAR, NMDAR traces. 

% move to correct start point
cd('E:\Experiments\NMDAR mEPSC\04_processing\Subtraction\2022 Analysis')

% choose whether to plot individual figures too
indi_fig = 0; % logical (1 plots, 0 skips)
if indi_fig == 1
    disp('Plotting individual figures')
else
    disp('Individual figures not plotted')
end

% run the loop three times as there's currently no E3M ... 
for j = 1:4 
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
    
    if indi_fig == 1
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
        scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor','yellow')
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
    else 
    end
    
    % return to home
    cd('E:\Experiments\NMDAR mEPSC\04_processing\Subtraction\2022 Analysis')

    % create output
    out(j).Condition = strrep(file(1:end-7),'_',' ');
    out(j).filepaths = paths;
    out(j).AMPAR = AMPAR;
    out(j).AMPAR_amp = AMPAR_amp;
    out(j).AMPAR_ave = AMPAR_ave;
    out(j).AMPAR_Hz = AMPAR_Hz;
    out(j).AMPAR_rise = AMPAR_rise;
    out(j).AMPAR_decay = AMPAR_decay;
    out(j).MIXED = MIXED;
    out(j).MIXED_amp = MIXED_amp;
    out(j).MIXED_ave = MIXED_ave;
    out(j).MIXED_Hz = MIXED_Hz;
    out(j).MIXED_rise = MIXED_rise;
    out(j).MIXED_decay = MIXED_decay;
    out(j).NMDAR = NMDAR;
    out(j).NMDAR_ave = NMDAR_ave;
    out(j).NMDAR_amp = NMDAR_amp;
    out(j).N_rec =  size(out(j).MIXED_amp,1);
    out(j).time = time;

   
    % clear variables
    clearvars -except out indi_fig time
end

%% create the mixed master figure
figure; 
    % Overlay (Traces)
    sgtitle('Mixed mEPSCs by Gender & Genotype')
    subplot(2,3,[1,4]); plot(time,out(1).MIXED_ave)
    hold on; plot(time,out(2).MIXED_ave); plot(time,out(3).MIXED_ave);
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)');
    xlim([0 0.1]);ylim([-25 5]); title('Ensemble Averages (Median) Overlaid Traces');
    legend(out(1).Condition,out(2).Condition,out(3).Condition,'location','southeast','linewidth',0.001);
    % Amplitude Summary
    subplot(2,3,2)
    Y = padcat(out(1).MIXED_amp,out(2).MIXED_amp,abs(out(3).MIXED_amp),abs(out(4).MIXED_amp));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition,out(4).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:4,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    scatter(x(:,4),Y(:,4),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0 0]);
    ylabel('Amplitude (pA)'); xlabel('Condition'); title('Amplitude')
    % Hz Summary
    subplot(2,3,3)
    Y = padcat(out(1).MIXED_Hz,out(2).MIXED_Hz,abs(out(3).MIXED_Hz),abs(out(4).MIXED_Hz));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition, out(4).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:4,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    scatter(x(:,4),Y(:,4),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0 0]);
    ylabel('Frequency (Hz)'); xlabel('Condition'); title('Frequency')
    % Rise Summary
    subplot(2,3,5)
    Y = padcat(out(1).MIXED_rise,out(2).MIXED_rise,abs(out(3).MIXED_rise), abs(out(4).MIXED_rise));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition, out(4).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:4,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    scatter(x(:,4),Y(:,4),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0 0]);
    ylabel('Rise Time Constant (Hz)'); xlabel('Condition'); title('Rise Time')
    % Decay Summary
    subplot(2,3,6)
    Y = padcat(out(1).MIXED_decay,out(2).MIXED_decay,abs(out(3).MIXED_decay),abs(out(4).MIXED_decay));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition, out(4).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:4,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    scatter(x(:,4),Y(:,4),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0 0]);
    ylabel('Decay Time Constant (Hz)'); xlabel('Condition'); title('Decay Time')
    
%% create the AMPAR master figure
figure; 
    % Overlay (Traces)
    sgtitle('AMPAR mEPSCs by Gender & Genotype')
    subplot(2,3,[1,4]); plot(time,out(1).AMPAR_ave)
    hold on; plot(time,out(2).AMPAR_ave); plot(time,out(3).AMPAR_ave); plot(time,out(4).AMPAR_ave);
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)');
    xlim([0 0.1]);ylim([-25 5]); title('Ensemble Averages (Median) Overlaid Traces');
    legend(out(1).Condition,out(2).Condition,out(3).Condition,'location','southeast','linewidth',0.001);
    % Amplitude Summary
    subplot(2,3,2)
    Y = padcat(out(1).AMPAR_amp,out(2).AMPAR_amp,abs(out(3).AMPAR_amp));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    ylabel('Amplitude (pA)'); xlabel('Condition'); title('Amplitude')
    % Hz Summary
    subplot(2,3,3)
    Y = padcat(out(1).AMPAR_Hz,out(2).AMPAR_Hz,abs(out(3).AMPAR_Hz));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    ylabel('Frequency (Hz)'); xlabel('Condition'); title('Frequency')
    % Rise Summary
    subplot(2,3,5)
    Y = padcat(out(1).AMPAR_rise,out(2).AMPAR_rise,abs(out(3).AMPAR_rise));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    ylabel('Rise Time Constant (Hz)'); xlabel('Condition'); title('Rise Time')
    % Decay Summary
    subplot(2,3,6)
    Y = padcat(out(1).AMPAR_decay,out(2).AMPAR_decay,abs(out(3).AMPAR_decay));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    ylabel('Decay Time Constant (Hz)'); xlabel('Condition'); title('Decay Time')
    
%% Create the NMDAR master figure
figure; 
    % Overlay (Traces)
    sgtitle('NMDAR mEPSCs by Gender & Genotype')
    subplot(2,3,[1,4]); plot(time,out(1).NMDAR_ave)
    hold on; plot(time,out(2).NMDAR_ave); plot(time,out(3).NMDAR_ave);
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (pA)'); xlabel('Time (s)');
    xlim([0 0.1]);ylim([-10 5]); title('Ensemble Averages (Median) Overlaid Traces');
    legend(out(1).Condition,out(2).Condition,out(3).Condition,'location','southeast','linewidth',0.001);
    % Amplitude Summary
    subplot(2,3,2)
    Y = padcat(abs(out(1).NMDAR_amp),abs(out(2).NMDAR_amp),abs(out(3).NMDAR_amp));
    b = boxplot(Y,...
        'labels',{out(1).Condition, out(2).Condition, out(3).Condition});
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    hold on
    x=repmat(1:3,length(Y),1);
    hold on
    scatter(x(:,1),Y(:,1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0 0.4470 0.7410]);
    scatter(x(:,2),Y(:,2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.8500 0.3250 0.0980]);
    scatter(x(:,3),Y(:,3),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',[0.9290 0.6940 0.1250]);
    ylabel('Amplitude (pA)'); xlabel('Condition'); title('Amplitude')
    for i = 1:3
        plot_array = [3 5 6];
        colour_array = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
        % Overlay loop
        subplot(2,3,plot_array(i)); plot(time,out(i).NMDAR,'color',[0.5 0.5 0.5 0.4],'HandleVisibility','off'); 
        hold on; plot(time,out(i).NMDAR_ave,'color', colour_array(i,:), 'linewidth' ,2)
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        ylabel('Amplitude (pA)'); xlabel('Time (s)');
        xlim([0 0.1]);ylim([-20 10]); title(out(i).Condition);
    end

    clearvars -except out