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
[file,path,~] = uigetfile('*.tdms');
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

%% Plots figures
% Check to see file loaded in correctly
if exist('S','var') == 1
    
    % Plot whole trace, filtered at 1kHz 
    Time = S.array(:,1);
    Trace = S.array(:,2);
    YF = filter1(Trace, Time, 0, 1000);
    figure;
    plot(YF)
    title('Filtered Whole Trace')
    hold on
    
    % data points 23606161 - 24421511 works nicely
    % y = 0 to -3 x10^-10 
    
    % Plot highlighted region on trace
    x1=23606161;
    x2=24421511;
    y1=0;
    y2=(-3*10^(-10));
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    plot(x, y, 'r-', 'LineWidth', 3);
    
    % Plot rolling trace
    figure;
    h = animatedline;
    numpoints = length(YF(23606161:24421511));
    axis([0 numpoints/sample_rate (-3*10^(-10)) 0])
    x = linspace(0,numpoints/sample_rate,numpoints);
    y = YF(23606161:24421511);
    a = tic; % start timer
    for k = 1:numpoints
        addpoints(h,x(k),y(k))
        b = toc(a); % check timer
        if b > (1/25)
            drawnow % update screen every 1/30 seconds
            a = tic; % reset timer after updating
        end
    end
    drawnow % draw final frame
    
end


 
%% Save recording as .phy
