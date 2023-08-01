%function SpikeAdapt

%% Decription
% Standaalone function to perform the spike frequency adaptation analysis
% as required in thesis corrections. To be implemented in IStep with
% successive updates.

%% Code

% close figures
close all

% determine the number of input arguments
%numInputs = nargin;

% load file of interest
disp('---------')
disp('Select the first file (from the first wave) in the recording, ephysIO will load the rest automatically')
disp('---------')
[file,path] = uigetfile('*.ma'); % select file of interest
clampfile = fullfile(path,file); % get full filepath of file
S = ephysIO(clampfile);
Time = S.array(:,1);
Waves = S.array(:,2:end);

% cd to path
splitPath = split(clampfile,filesep);
newPath = char(string(join(splitPath(1:end-2),filesep)));
cd(newPath)
%% run analysis or not?
% gives the user the option to abort if the data looks horrendous

% plot overall figure
fh = figure();
plot(Time,Waves(:,1)*1000,'color','black')
hold on; plot(Time,Waves*1000,'color','black','HandleVisibility','off')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
ylabel('Membrane Potential (mV)'); xlabel('Time (s)')
title('Current Step Waveform')

% ask the user 
dlgTitle    = 'Run Analysis';
dlgQuestion = 'Would you like to run the Current Step Analysis?';
run = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');

if run == "Yes"
    disp('---------')
    disp('Performing Analysis, please wait ...')
    disp('---------')

    % select region to search for action potentials
    title('Select detection region')
    [detection,~] = ginput(2);
    detStart = round(detection(1)/S.xdiff); % start of detection in data points
    detEnd = round(detection(2)/S.xdiff); % end of detection in data points
    xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off')
    xline(detection(2),'linestyle','--','color','blue','linewidth',2)

    % select region to detect the baseline from
    title('Select baseline region')
    [baseline,~] = ginput(2);
    baseStart = round(baseline(1)/S.xdiff); % start of detection in data points
    baseEnd = round(baseline(2)/S.xdiff); % end of detection in data points
    xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off')
    xline(baseline(2),'linestyle','--','color','green','linewidth',2)
    legend('Recording Data','AP Detection Region','Baseline Detection Region','linewidth',1)
    title('Detection Regions')
    Vbase = mean(Waves(baseStart:baseEnd,:)); % determine the baseline (intra-step)

    % pause for 2 seconds to allow user to visualise regions, then close
    pause(2)
    % close figure
    close(fh)

        % preallocate pks, locs, w, p, numSpikes and mute error here
    warning('off','signal:findpeaks:largeMinPeakHeight');
    pks = cell(size(Waves,2),1);
    locs = cell(size(Waves,2),1);
    difflocs = cell(size(Waves,2),1); % difference between loc in dp
    normlocs = cell(size(Waves,2),1); % difference between loc in dp
    w = cell(size(Waves,2),1);
    p = cell(size(Waves,2),1);
    numSpikes = zeros(size(Waves,2),1);
    nlocs = cell(size(Waves,2),1); % interval index
    
    
    % findpeaks to determine number and location of AP's per wave
    for i = 1:size(Waves,2)
        % pks = value of peak
        % locs = index of peak
        % w = width
        % p = prominence
        % hard coded to only look for AP's between detStart and detEnd
        [pks{i},locs{i},w{i},p{i}] = findpeaks(Waves(round(detStart):round(detEnd),i),...
            'MinPeakHeight',0,...
            'MinPeakProminence',0.01,...
            'MaxPeakWidth',0.02*10^4,...
            'MinPeakDistance',0.02*10^4,...
            'WidthReference','halfheight');
        numSpikes(i) = size(pks{i},1);
        difflocs{i} = diff(locs{i}); % difference between loc in dp
        normlocs{i} = normalize(difflocs{i},'scale','first'); % difflocs norm to first value
        nlocs{i} = 1:size(difflocs{i}); % interval index
    end
    
    % plot whole of current step protocol
    fh = figure();
    fh.WindowState = 'maximized'; subplot(7,4,[1,5]); plot(Time,Waves*1000,'color','black','HandleVisibility','off')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Membrane Potential (mV)');
    title('Current Step Waveform')
    ax = gca; xax = ax.XAxis; set(xax,'visible','off')
    hold on
    xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
    xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
    xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
    xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
    hold off
    legend('Detection Region','Baseline Region','linewidth',1)

     % plot waveform of current step read from second channel of clampfit
    SI = ephysIO({clampfile,2}); % load in the current data
    x = SI.array(:,1); % x is time here
    pA_waveform = SI.array(:,2:end); % y is the array of current data
    Ibase = mean (mean(pA_waveform(baseStart:baseEnd,:))); % determine the baseline (intra-step)
    N = size(SI.array,2) - 1; % get the number of waves in the array (minus time)
    lo = dsearchn(x,detection(1)) + 1; % start of the test pulse
    hi = dsearchn(x,detection(2)) - 1; % end of the test pulse
    pA = zeros(1,N); % preallocate a blank series of steps
    for i = 1:N
       pA(i) = mean (pA_waveform(round(lo+(detEnd-detStart)*0.1):...
           round(hi-(detEnd-detStart)*0.1),i)); % fill the steps with the mean of each step
    end
    pA = fix((pA - Ibase) * 1e+12); % round to zero, baseline subtract and put into pA
    subplot(7,4,9); plot(x,(pA_waveform-Ibase)*1e12,'linewidth',1,'color','black','HandleVisibility','off')
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    xlabel('Time (s)'); ylabel('Command(pA)'); ylim([min(min(pA_waveform*1e12))-50,max(max(pA_waveform*1e12))+50])
    ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')
    hold on
    xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
    xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
    xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
    xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
    hold off
    
    % First AP of interest values
    C = locs;
    idx = ~cellfun('isempty',C);
    outs = zeros(size(C));
    outs(idx) = cellfun(@(v)v(1),C(idx));
    for lp = 1:size(outs,1)
        if outs(lp) > 0
            outs(lp) = outs(lp) + detStart; % account for the detection zone
        end
    end
    subplot(7,4,[2,6,10]);
    plot(Time,Waves(:,end-6)*1000,'color','red'); hold on; plot(Time,Waves(:,11)*1000,'color','black')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
    lgd = legend(char(string(pA(end-6))),char(string(pA(11))),'linewidth',1);
    title(lgd,'Current (pA)')
    title('Exemplary Waves')

    % plot locations of peaks
    
else 
end