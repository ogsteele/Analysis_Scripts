%% NSFA Selector

% What do we need the script to do?
% logic should be ...
    % plot trace 
    % ginput to square off the region you want the peaks to be within
    % series of logical statements
        % does this wave have a peak in that square region?
        % if yes, merge into new recording
        % if no, discard
    % save the new recording for use in the NSFA script

        
% load the file with ephysIO and plot
quickfig

% close the second fig that opens
close(gcf)

% adjust axis label to make sense
ylabel('Current (A)')

% close spare fig and zoom into the area of interest and hit enter
disp('Close spare figure if needed and zoom into trace, before hitting Enter')
pause % pause until key press

% select the baseline region for baseline subtraction
disp('Select the baseline region for baseline subtraction')
[base_s,~] = ginput(2);
xline(base_s(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off')
xline(base_s(2),'linestyle','--','color','green','linewidth',2,'DisplayName','Baseline Region')
legend()

% calculate baseline values
base_x = round(base_s/S1.xdiff);
baseVal = mean(Waves(base_x(1):base_x(2),:));

% close spare fig and zoom into the area of interest and hit enter
disp('Close spare figure if needed and zoom into trace, before hitting Enter')
pause % pause until key press

% select the total spike region
disp('Select the total spike region')
[spike_x,~] = ginput(2);
xline(spike_x(1),'linestyle','--','color','red','linewidth',2,'HandleVisibility','off')
xline(spike_x(2),'linestyle','--','color','red','linewidth',2,'DisplayName','Detection Region')
legend()

% ginput to select points and create square
disp('Zoom in further if needed')
pause

disp('Select point one side of peak of interest')
[x1,y1] = ginput(1); % select upper left of peak of interest
h(1) = xline(x1,'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off');
h(2) = yline(y1,'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off');
disp('Select opposite point to close the square')
[x2,y2] = ginput(1); % select lower right of peak of interest
h(3) = xline(x2,'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off');
h(4) = yline(y2,'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off');

% plot square of selection zone
sq_x = [x1, x2, x2, x1, x1];
sq_y = [y1, y1, y2, y2, y1];
hold on; plot(sq_x, sq_y,...
    'linestyle','--',...
    'color','blue',...
    'linewidth',2,...
    'DisplayName','Selection Region');

delete(h(1:4))

% organise points into square
select_x = [x1 x2]; y = [y1 y2]; % concatenate arrays for below admin
northWest = [min(select_x) max(y)];
northEast = [max(select_x) max(y)];
southWest = [min(select_x) min(y)];
southEast = [max(select_x) min(y)];

% convert x and spikes into datapoints
select_x = floor(select_x/S1.xdiff);
spike_x = floor(spike_x/S1.xdiff);

% preallocate the logical loop
peakLogical = zeros(size(Waves,2),1);
ind = zeros(size(Waves,2),1);
val = zeros(size(Waves,2),1);

% logical tests in for loop per wave
for i =  1:size(Waves,2)
    % ind ~= x values, val ~= y
    [val(i),ind(i)] = min(Waves(min(spike_x):max(spike_x),i));
    % ind(i) > min(x)
    % ind(i) < max(x)
    % val(i) > min(y)
    % val(i) < max(y)
    ind(i) = ind(i) + (min(spike_x)-2); % correct for window offset
    if ind(i) > min(select_x) && ind(i) < max(select_x) ... 
            && ...
       val(i) > min(y) && val(i) < max(y)
        peakLogical(i) = true;
        hold on; plot(ind(i)*S1.xdiff,val(i),'*r','HandleVisibility','off')
    else
        peakLogical(i) = false;
    end
end
hold on; pc = plot(ind(1)*S1.xdiff,val(1),'*r','DisplayName','Peak Amplitude');

% convert peakLogical to a logical
peakLogical = logical(peakLogical)';

% extract only the waves that pass the double logic gates above
selectWaves = Waves(:,peakLogical);

% pull xlim and ylim from gcf
yl = ylim; xl = xlim;
title('NSFA Selector.m output');

% replot with old and new in different styles
delete(fw)
plot(Time,Waves,'color','black','linestyle','--','HandleVisibility','off')
plot(Time,selectWaves,'color','black','linewidth',1.5,'HandleVisibility','off')
plot(nan, nan, 'color', 'black', 'linestyle', '--','DisplayName','Deleted Traces');
plot(nan, nan, 'color', 'black', 'linewidth', 2, 'DisplayName','Selected Traces');

% perform baseline subtraction here 
% (note, it won't appear in current fig but does affect the math and the
% final trace produced and saved)
for i = 1:size(Waves,2)
    Waves(:,i) = Waves(:,i) - baseVal(i);
end

% replot with the circles on the top
delete(pc)
for i = 1:size(peakLogical,2)
    if peakLogical(i) == true
        plot(ind(i)*S1.xdiff,val(i),'*r','linewidth',2,'HandleVisibility','off')
    end
end
L1 = find(peakLogical, 1, 'first');
plot(ind(L1)*S1.xdiff,val(L1),'*r','linewidth',2,'DisplayName','Peak Amplitude')

% create variables
wavenum = 1:size(val,1);
rawPeaks = val*10e9;
selectPeaks = val(peakLogical)*10e9;
selectWavesind = find(peakLogical);

% create the linear fit
x = selectWavesind;
y = selectPeaks;
p = polyfit(x,y,1);
f = polyval(p,x);

% plot variables now
figure; 
plot(wavenum,rawPeaks,'o','color',[0.5 0.5 0.5],'DisplayName','deleted peaks')
hold on; h(5) = plot(selectWavesind,selectPeaks,'o','color','red','LineWidth',2,'DisplayName','selected peaks');
hold on; plot(selectWavesind,movmean(selectPeaks,4),'DisplayName','movmean (n = 4)')
hold on; plot(x,f,'-','DisplayName','Linear Fit')

% format the plot
box off; set(gcf,'color','white'); set(gca,'linewidth',2);
xlabel('Wave number'); ylabel('Current (nA)'); title('Amplitude as a function of time');
legend()

% zoom in if needeed
disp('zoom into selected waves if needed, before hitting Enter')
pause % pause until key press

% select the total spike region
disp('Select the region w/o systematic change')
[selectRegion,~] = ginput(2);
xline(selectRegion(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off')
xline(selectRegion(2),'linestyle','--','color','blue','linewidth',2,'DisplayName','Selected Region')
legend()

% select the waves you now want to keep
waveLogical = wavenum((wavenum > selectRegion(1)) ...
                        & (wavenum < selectRegion(2))); % above below lines
waveLogical = waveLogical(...
                    peakLogical(waveLogical(1):waveLogical(end))); % is peak logical?
selectPeaks = val(waveLogical)*10e9;

% update the plot
delete(h(5))
hold on; plot(waveLogical,selectPeaks,'o','color','red','LineWidth',2,'DisplayName','selected peaks');

% extract only the waves that pass the double logic gates above
selectWaves = Waves(:,waveLogical);

% save the new file (with full t)
cd(path); cd ..; % navigate to one dir above the file recs
splitfile = split(clampfile,filesep);
savefile = [char(splitfile(end-2)) '.phy'];
new_array = [Time,selectWaves];
ephysIO(savefile,new_array,S1.xunit,S1.yunit)

% save another file that only has the region of interest
savefile = [char(splitfile(end-2)) '_ROI.phy'];
start_t = round(0.178/S1.xdiff);
end_t = round(0.195/S1.xdiff); 
w = selectWaves(start_t:end_t,:);
t = Time(1:size(w,1));
ephysIO(savefile,[t w],S1.xunit,S1.yunit)

% save the first figure
saveFigName = char(splitfile(end-2));
savefig(saveFigName)

% plot the ROI trace
figure; plot(t,w,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2);
xlabel('Time (s)'); ylabel('Current (A)'); title('Selected Waves');

% print out savefile so it's the last thing in the command window
disp(savefile)