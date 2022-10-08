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
xlabel('Current (A)')

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
hold on; plot(ind(1)*S1.xdiff,val(1),'*r','DisplayName','Peak Amplitude')

% convert peakLogical to a logical
peakLogical = logical(peakLogical)';

% extract only the waves that pass the logic gate above
selectWaves = Waves(:,peakLogical);

% pull xlim and ylim from gcf
yl = ylim; xl = xlim;
title('Raw Waves');

% plot select waves into a new graph on the same axis limits
figure; 
plot(Time,selectWaves,'color','black')
xlim(xl); ylim(yl);
box off; set(gcf,'color','white'); set(gca,'linewidth',2);
xlabel('Time (s)'); ylabel('Current (A)'); title('Selected Waves');

% save the new file
cd(path); cd ..; % navigate to one dir above the file recs
splitfile = split(clampfile,filesep);
savefile = [char(splitfile(end-2)) '.phy'];
new_array = [Time,selectWaves];
ephysIO(savefile,new_array,S1.xunit,S1.yunit)
