%% AutoPlot script
% Script to plot the data on a rolling window in 'real time' as determined
% by the sampling rate of the recording

%% Parameters
sample_rate = 20000;

%% Load in Data

% Select raw trace to visualise
title_str = "1. Select raw file of recording";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Select raw .tdms of recording');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    % Navigate to directory and load file with ephysIO
    cd(path)
    S = ephysIO(file);
end
clear('file', 'path')
%% Plots figures
% Check to see file loaded in correctly
if exist('S','var') == 1
    close all
    % Plot whole trace, filtered at 1kHz 
    Time = S.array(:,1); 
    Trace = S.array(:,2); 
    x_Time = linspace(0,length(Time)/sample_rate,length(Time)); % x axis
    YF = filter1(Trace, Time, 0, 500); % filter data in Hz
    YF_pA = YF*10^12; % scale data to pA
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    plot(x_Time,YF_pA)
    
    % style plot
    title('Filtered Whole Trace','FontSize',20)
    ylabel('Current (pA)')
    xlabel('Time (s)')
    set(gcf,'color','w');
    box off
    set(gca,'linewidth',4)
    
    % data points 23606161 - 24421511 works nicely
    % y = 0 to -3 x10^-10 
    
    % Plot rolling trace
    subplot(2,1,2)
    % style plot
    set(gcf,'color','w');
    box off
    set(gca,'linewidth',4)
    title('"Live" view of mEPSC acqusition','FontSize',20)
    ylabel('Current (pA)')
    xlabel('Time (s)')   
    h = animatedline;
    numpoints = 20*sample_rate;
    axis([0 numpoints/sample_rate -60 0]) % initially set to empty
    
    % get points from the first graph
    subplot(2,1,1)
    %[xi,~] = ginput(1);
    [xi] = getpts;
    selected_point = round(xi*sample_rate);
    
    % Plot highlighted region on trace
    hold on
    subplot(2,1,1)
    x1=selected_point/sample_rate; % convert point to seconds
    x2=selected_point/sample_rate+20;
    y1=+60;
    y2=(-60);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    plot(x, y, 'r-', 'LineWidth', 3);
    hold off
    
    hold on
    subplot(2,1,2)
    axis([0 numpoints/sample_rate ...
        (median(YF_pA(selected_point:selected_point+10000))-60) ...
        (median(YF_pA(selected_point:selected_point+10000))+60)]) % change to reflect point
    
    x = linspace(0,numpoints/sample_rate,numpoints);
    y = YF_pA(selected_point:selected_point+(20*sample_rate));
    a = tic; % start timer
    for k = 1:numpoints
        addpoints(h,x(k),y(k))
        b = toc(a); % check timer
        if b > (1/15) % update screen every 1/30 seconds
            drawnow 
            a = tic; % reset timer after updating
        end
    end
    drawnow % draw final frame
    
    %h = animatedline;
    %numpoints = length(YF_pA(23606161:24421511));
    %axis([0 numpoints/sample_rate -60 0])
    %x = linspace(0,numpoints/sample_rate,numpoints);
    %y = YF_pA(23606161:24421511);
    %a = tic; % start timer
    %for k = 1:numpoints
    %    addpoints(h,x(k),y(k))
    %    b = toc(a); % check timer
    %    if b > (1/15) % update screen every 1/30 seconds
    %        drawnow 
    %        a = tic; % reset timer after updating
    %    end
    %end
    %drawnow % draw final frame
end

% Tidy up
clear('h', 'numpoints','x','y','a','k','b',...
    'x1','y1','x2','y2','YF','YF_pA',...
    'Trace','Time','x_Time',...
    'amplifier_gain','amplifier_scale','ans','S','sample_rate')   

function [tprm,med,x,x2] = TestPulse_rm(data,parameters)

% remove test pulse and plot without
% calculate the pre pulse mean to fill the gaps
tprm_splits = data;
base = zeros(size(tprm_splits,2),1);
median_YF = zeros(size(tprm_splits,2),1);
for i = 1:size(tprm_splits,2)
    base(i) = trimmean(tprm_splits(1:parameters.pulse_start,i),33,'floor'); % baseline average (in pA)
    tprm_splits(parameters.pulse_start:(parameters.pulse_end+250),i) = base(i); % swap pulse for base
    median_YF(i) = median(tprm_splits(:,i)); % calculate the median of the splits
end
% concatenate tprm_splits
tprm_splits_conc = vertcat(tprm_splits(:));
% apply a filter to clean the look of the data
t = (0:size(tprm_splits_conc,1)-1)';
t = t./parameters.sample_rate;
YF = filter1(tprm_splits_conc, t, 0, 300);
x = linspace(0,(size(tprm_splits_conc,1)/20000),size(tprm_splits_conc,1));
%figure
%plot(x,YF,'Color',[0 0.4470 0.7410 0.2]) % to make lighter
%title('filter')
x2 = linspace(0,(size(tprm_splits_conc,1)/20000),size(median_YF,1));
%hold on
%plot(x2,smooth(median_YF),'linewidth',4,'color',[0 0.4470 0.7410]) % not made lighter
tprm = YF;
med = median_YF;
end
