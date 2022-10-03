% devThresh
% development version of threshold detector taken from Istep (26.09.22)

% load file of interest
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

% select region to search for action potentials
detStart = 0.5/S.xdiff; % start of detection in data points
detEnd = 1.5/S.xdiff; % end of detection in data points

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
logicalIndexes =  outs < 27500 & outs > 1;
wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on

%% AP analysis - Surpasses 4 mV/ms

% determine window size
% 25 is good for organotypic
% 50 is good for slower APs
winSize_ms = inputdlg('Choose AP window size (ms)','winSize_ms',1,{'longer window for slow APs'});
winSize_ms = str2num(cell2mat(winSize_ms)); % convert to number
winSize_ms = (winSize_ms*1e-3)/S.xdiff; % (ms to s)to data points


% plot action potential
warning('off','MATLAB:colon:nonIntegerIndex')
AP_Window = Waves(outs(wavenum_first)-(winSize_ms/2):outs(wavenum_first)+(winSize_ms/2),wavenum_first)*1000;
warning('on','MATLAB:colon:nonIntegerIndex')
% % 
% % % Overshoot in mV
% % [Overshoot,ind_o] = max(AP_Window);
% % % Afterhyperpolarisation in mV
% % [Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% % % Baseline
% % Base = mean(AP_Window(1:350));
% % 
% % % Action potential halfwidth
% % % Halfwidth in ms
% % figure;
% % % subplot(7,4,[4,8,12]);
% % [~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
% %         'MaxPeakWidth',0.02*10^4,...
% %         'MinPeakDistance',600,...
% %         'WidthReference','halfheight',...
% %         'Annotate','extent');
% % findpeaks(AP_Window-Base,'MinPeakHeight',0,...
% %         'MaxPeakWidth',0.02*10^4,...
% %         'MinPeakDistance',600,...
% %         'WidthReference','halfheight',...
% %         'Annotate','extent');
% % box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
% % set(gca,'linewidth',2); set(gcf,'color','white'); title('4 mV/ms');
% % ylim([-40 120])
% % Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms
% % 
% % % Action Potential Threshold
% % 
% % % dv/dt surpasses 4mV/ms (dv/dt > 1)
% % N = diff(AP_Window(1:ind_o-20)); % period start to just before the peak
% % [closestValue, closestIndex] = min(abs(N - 1.'));
% % hold on; plot((diff(AP_Window))*5-20)
% % hold on; plot(closestIndex,AP_Window(closestIndex)-Base,'or')
% % ind_t = closestIndex; % threshold index
% % Threshold = AP_Window(ind_t); % in mV
% % 
% % hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after
% % 
% % % recalculate (and overwrite) the amplitude now you have a threshold
% % Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV
% % 
% % % Depolarisation Rate
% % % identify the closest value (and index) to the threshold value after peak
% % N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
% % [~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
% % [~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));
% % 
% % % between 20 % and 80 % of the rise phase (based on amplitude, not index)
% % rise_20_ind = closestIndex_20 + ind_t;
% % rise_80_ind = closestIndex_80 + ind_t;
% % Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
% % hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)
% % 
% % % Repolarisation Rate
% % 
% % % identify the closest value (and index) to the threshold value after peak
% % N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
% % [~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
% % [~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));
% % 
% % % between 20 % and 80 % of the rise phase (based on amplitude, not index)
% % fall_20_ind = closestIndex_80 + ind_o;
% % fall_80_ind = closestIndex_20 + ind_o;
% % Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
% % hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)
% % 
% % legend('trace','peak', ...
% %     'amplitude','halfwidth', ...
% %     'border','threshold', ...
% %     'dv/dt',...
% %     'hyperpolar.',...
% %     'rise','fall',...
% %     'Location','northeast')
% % 
% % savefig('method1.fig')
% % pause(2)
% % %% AP analysis - local min of the peak of the third derivative
% % 
% % 
% % % Overshoot in mV
% % [Overshoot,ind_o] = max(AP_Window);
% % % Afterhyperpolarisation in mV
% % [Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% % % Baseline
% % Base = mean(AP_Window(1:350));
% % 
% % % Action potential halfwidth
% % % Halfwidth in ms
% % figure;
% % % subplot(7,4,[4,8,12]);
% % [~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
% %         'MaxPeakWidth',0.02*10^4,...
% %         'MinPeakDistance',600,...
% %         'WidthReference','halfheight',...
% %         'Annotate','extent');
% % findpeaks(AP_Window-Base,'MinPeakHeight',0,...
% %         'MaxPeakWidth',0.02*10^4,...
% %         'MinPeakDistance',600,...
% %         'WidthReference','halfheight',...
% %         'Annotate','extent');
% % box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
% % set(gca,'linewidth',2); set(gcf,'color','white'); title('local min, total max, triple ndiff');
% % ylim([-40 120])
% % Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms
% % 
% % % Action Potential Threshold
% % 
% % % triple ndiff
% % diffWin = AP_Window;
% % for diffnum = 1:3
% %     diffWin = ndiff(diffWin,Time(1:size(diffWin,1)));
% % end
% % 
% % % plot the differntial below the trace
% % y = smooth(diffWin)*5e-12-20;
% % hold on; plot(y); % plotting of the actual triple diff trace
% % x = 1:size(y,1);
% % TF = islocalmin(y); % works for the slower action potential, not faster
% % % perhaps consider 'MinSeperation' val as a function of speed of action
% % % potential? Maybe related halfwidth.
% % 
% % % discover the local min to the left of the total max
% % TF_ind = find(TF); % return array of all the local minima index
% % [val,ind] = max(y(1:end-5)); % find the max of y
% % neg_idx = (TF_ind - ind) < 0; % return a logical index of the the negative (left shifted) min indexes
% % max_neg = max(TF_ind(neg_idx)); % return the largest index to the left of the overall maximum
% % ind_t = max_neg; % threshold index
% % Threshold = AP_Window(ind_t); % in mV <-- to go to the output
% % hold on; plot(ind_t,y(ind_t),'*b') % plot the local min
% % hold on; plot(ind_t,AP_Window(ind_t)-Base,'or') % plot the threshold
% % 
% % % plot the afterhyperpolarisation
% % hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after
% % 
% % % recalculate (and overwrite) the amplitude now you have a threshold
% % Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV
% % 
% % % Depolarisation Rate
% % % identify the closest value (and index) to the threshold value after peak
% % N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
% % [~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
% % [~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));
% % 
% % % between 20 % and 80 % of the rise phase (based on amplitude, not index)
% % rise_20_ind = closestIndex_20 + ind_t;
% % rise_80_ind = closestIndex_80 + ind_t;
% % Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
% % hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)
% % 
% % % Repolarisation Rate
% % 
% % % identify the closest value (and index) to the threshold value after peak
% % N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
% % [~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
% % [~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));
% % 
% % % between 20 % and 80 % of the rise phase (based on amplitude, not index)
% % fall_20_ind = closestIndex_80 + ind_o;
% % fall_80_ind = closestIndex_20 + ind_o;
% % Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
% % hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)
% % 
% % legend('trace','peak', ...
% %     'amplitude','halfwidth', ...
% %     'border', ...
% %     'triple ndiff','local minima','threshold',...
% %     'hyperpolar.',...
% %     'rise','fall',...
% %     'Location','northeast')
% % 
% % savefig('method2.fig')
% % pause(2)
%% AP analysis - initial peak of the second derivative


% Overshoot in mV
[Overshoot,ind_o] = max(AP_Window);
% Afterhyperpolarisation in mV
[Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% Baseline
Base = mean(AP_Window(1:350));

% Action potential halfwidth
% Halfwidth in ms
figure;
% subplot(7,4,[4,8,12]);
[~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
set(gca,'linewidth',2); set(gcf,'color','white'); title('first peak of double ndiff');
ylim([-40 120])
Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms

% Action Potential Threshold

% triple ndiff
diffWin = AP_Window;
for diffnum = 1:2
    diffWin = ndiff(diffWin,Time(1:size(diffWin,1)));
end

% LPF_Hz Note <-- faster the action potential, lower the threshold
% 330 works well for organotypic, 1000 works well for iPSCs

LPF_Hz = inputdlg('Choose LPF cut off (Hz)','LPF_Hz',1,{'lower cutoff for fast APs'});
LPF_Hz = str2num(cell2mat(LPF_Hz));

% plot the differntial below the trace
y = filter1(diffWin,Time(1:size(diffWin,1)),0,LPF_Hz)*5e-8-20; 
hold on; plot(y); % plotting of the actual triple diff trace
x = 1:size(y,1);

% find peaks in the double diff trace
[dd_pks,dd_locs,~,~] = findpeaks(y,"MinPeakProminence",0.5); 
% make sure to understand peak prominence properly

% discover the local min to the left of the total max
ind_t = dd_locs(1); % threshold index is the first peak detected
Threshold = AP_Window(ind_t); % in mV <-- to go to the output
hold on; plot(ind_t,y(ind_t),'*b') % plot the initial peak
hold on; plot(ind_t,AP_Window(ind_t)-Base,'or') % plot the threshold

% plot the afterhyperpolarisation
hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after

% recalculate (and overwrite) the amplitude now you have a threshold
Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV

% Depolarisation Rate
% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
rise_20_ind = closestIndex_20 + ind_t;
rise_80_ind = closestIndex_80 + ind_t;
Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)

% Repolarisation Rate

% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
fall_20_ind = closestIndex_80 + ind_o;
fall_80_ind = closestIndex_20 + ind_o;
Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)

legend('trace','peak', ...
    'amplitude','halfwidth', ...
    'border', ...
    'double ndiff','initial max','threshold',...
    'hyperpolar.',...
    'rise','fall',...
    'Location','northeast')

savefig('method3.fig')