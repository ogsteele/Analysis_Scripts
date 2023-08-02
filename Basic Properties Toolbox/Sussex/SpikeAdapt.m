function [output] = SpikeAdapt(I,auto,path)
%% Decription
% Standalone function to perform the spike frequency adaptation analysis
% as required in thesis corrections. To be implemented in IStep with
% successive updates. Auto and path behaviour are more specific to speed
% analysing and creating loops, redudant if not needed.

% SpikeAdapt v1 (last updated: 02/08/23)
% Author: OGSteele

% example use;
%       output = SpikeAdapt(I,path);
%       for automated running
%       eg. adaptation = SpikeAdapt(300,1,'E:\Experiments\Ketamine mEPSC\04_processing\20211211\Slice 1\Cell 2\CC_IV_1_full_001\000\Clamp1.ma')

% output 
%   is a matlab structure containing the following;
%       filepath - filepath of the first wave
%       peakTimes - time of peaks in s
%       peakAmps - overshoot of peaks in mV
%       ISI - interspike interval in s
%       peakHz - instantaneous spike frequency in Hz
%       adaptIndex - adaptation index

% input
%       I is the current value at which to measure the spike accomodation
%           in pA (eg. 300)
%       auto (1 or 0) is a logical decision of whether or not to
%           automatically select the recording or to manually select. If
%           manual (auto = 0), UI input arises. If auto (auto = 1), also specify a path
%       path is the filepath of the first clampfile

%% Code


% close figures
close all

% determine the number of input arguments
numInputs = nargin;

if numInputs == 0 
    % if the numput of input arguments is zero, assume manual selection
    % with default parameters
    I = 300;
    auto = 0;
    path = '';
end

% determine whether automated selection or manual selection
if auto == 1
    clampfile = path;
    S = ephysIO(clampfile);
    Time = S.array(:,1);
    Waves = S.array(:,2:end);
    % cd to path
    splitPath = split(clampfile,filesep);
    newPath = char(string(join(splitPath(1:end-2),filesep)));
    cd(newPath)
    run = "Yes";
    fh = figure();
    plot(Time,Waves(:,1)*1000,'color','black')
    hold on; plot(Time,Waves*1000,'color','black','HandleVisibility','off')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Membrane Potential (mV)'); xlabel('Time (s)')
    title('Current Step Waveform')
elseif auto == 0

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
    
end



if run == "Yes"
    disp('---------')
    disp('Performing Analysis, please wait ...')
    disp('---------')

    % select region to search for action potentials
    title('Select detection region')
    detection = [0.4917;1.5230];
    detStart = round(detection(1)/S.xdiff); % start of detection in data points
    detEnd = round(detection(2)/S.xdiff); % end of detection in data points
    xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off')
    xline(detection(2),'linestyle','--','color','blue','linewidth',2)

    % select region to detect the baseline from
    title('Select baseline region')
    baseline = [0.1381;0.4401];
    baseStart = round(baseline(1)/S.xdiff); % start of detection in data points
    baseEnd = round(baseline(2)/S.xdiff); % end of detection in data points
    xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off')
    xline(baseline(2),'linestyle','--','color','green','linewidth',2)
    legend('Recording Data','AP Detection Region','Baseline Detection Region','linewidth',1)
    title('Detection Regions')
    Vbase = mean(Waves(baseStart:baseEnd,:)); % determine the baseline (intra-step)

    % pause for 2 seconds to allow user to visualise regions, then close
    pause(1)
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
    fh.WindowState = 'maximized'; subplot(7,3,[1,4]); plot(Time,Waves*1000,'color','black','HandleVisibility','off')
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
    subplot(7,3,7); plot(x,(pA_waveform-Ibase)*1e12,'linewidth',1,'color','black','HandleVisibility','off')
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    xlabel('Time (s)'); ylabel('Command(pA)'); ylim([min(min(pA_waveform*1e12))-50,max(max(pA_waveform*1e12))+50])
    ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')
    hold on
    xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
    xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
    xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
    xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
    hold off

    % catch for scaling error
    if sum(pA) > 10000
        pA = [-200:20:400]
    else
    end

    % Plot the wave of interest, in our case closest to 300 pA stimulation
    stimVal = I; % change value if required
    [val,ind] = min(abs(pA-stimVal));
    subplot(7,3,[2,5,8]);
    plot(Time,Waves(:,ind)*1000,'color','red'); hold on; plot(Time,Waves(:,11)*1000,'color','black')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
    lgd = legend(char(string(pA(ind))),char(string(pA(11))),'linewidth',1);
    title(lgd,'Current (pA)')
    title('Exemplary Waves')

    % extract and plot peak coordinates
    peakInd = cell2mat(locs(ind));
    peakVal = cell2mat(pks(ind));
    peakTimes = zeros(size(peakInd));
    peakAmps = zeros(size(peakVal));
    for pI = 1:size(peakTimes,1)
        peakTimes(pI) = Time(peakInd(pI)+detStart);
        peakAmps(pI) = peakVal(pI)*1000;
    end
    plot(peakTimes,peakAmps,'--ob','DisplayName','Peak Intervals')

    % convert to instantaneous frequency and plot
    ISI = diff(peakTimes);
    peakHz = 1./ISI;
    subplot(7,3,[3,6,9])
    plot(peakHz,'-o',...
        'LineWidth',1,...
        'MarkerSize',10,...
        'MarkerEdgeColor','r')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    xlabel('Interspike Interval #'); ylabel('Instantaneous Frequency (Hz)');
    legend('Spike Frequency','linewidth',1)

    % calculate adaptation index and plot
    adaptIndex = ISI(1)/ISI(end);
    subplot(7,3,[13,16,19]);
    names = char('','Initial','Final','');
    x = [2,3];
    y = [ISI(1),ISI(end)];
    plot(x,y,'--ob',...
        'LineWidth',1.5,...
        'MarkerSize',6,...
        'MarkerEdgeColor','r')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    set(gca,'XTick',1:4,'XTickLabel',names)
    xlim([1,4]);
    xlabel('Interspike Interval #'); ylabel('Insterspike Interval (s)');
    txt = {['Adaptation Index: '],[num2str(adaptIndex)]};
    annotation('textbox','String',txt,'Position',subplot(7,3,[13,16,19]).Position,'VerticalAlignment','top','HorizontalAlignment','right','FitBoxToText','on','linewidth',1.5);

    % plot as inst freq over time with raster and current injection
    subplot(7,3,[14:15,17:18,20:21]);
    fTime(1:2) = [0,peakTimes(1)-0.0001];
    for fT = 1:size(peakTimes,1)
        fTime(2+fT) = peakTimes(fT);
    end
    fVals(1:3) = [0,0,0];
    for fV = 1:size(peakHz,1)
        fVals(3+fV) = peakHz(fV);
    end
    plot(fTime,fVals,'linewidth',1.5,'DisplayName','Instantaneous Frequency (Hz)')
    xlabel('Time (s)'); ylabel('Instantaneous Frequency (Hz)');
    hold on; plot(peakTimes(1),65,"|",'color','r','linewidth',3,'MarkerSize',10,'DisplayName','Spikes')
    hold on; plot(peakTimes,65,"|",'color','r','linewidth',3,'MarkerSize',10,'HandleVisibility','off')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    x = [0.5,1.5]; y = [-10,-10]; hold on; line(x,y,'linewidth',10,'color','black','DisplayName','Current Injection')
    legend('linewidth',1)
    xlim([0.4,1.4]); ylim([-20,70])

    %% create output and save
    output.filepath = path;
    output.peakTimes = peakTimes;
    output.peakAmps = peakAmps;
    output.ISI = ISI;
    output.peakHz = peakHz;
    output.adaptIndex = adaptIndex;

    f = waitbar(0,'Saving adapt figure ...');
    set(f,'Name','Saving output, do not close figures');
    outname = split(strtrim(clampfile),filesep);
    outname = char(string(outname(end-2))); 
    % navigate to root dir
    cd(newPath) 
    saveas(fh,[outname,'_adapt.fig']); % save the master fig
        waitbar(.5,f,'Saving data structure ...');
    outname = [outname '_adapt.dat']; % name the output the same as the folder the recording came from
    save([outname,'.mat'],'output')
        waitbar(1,f,'Finishing');
        close(f)
else 
end