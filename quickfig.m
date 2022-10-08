% quickfig script to quickly plot and visualise recordings w/ ephysIO

% select the file of interest and load with ephysIO
[file,path] = uigetfile();
clampfile = fullfile(path,file);
S1 = ephysIO({clampfile,1}); % channel one
S2 = ephysIO({clampfile,2}); % channel two

% split into waves of interest
Time = S1.array(:,1);
Waves = S1.array(:,2:end);
Command = S2.array(:,2:end);

% plot figures
figure; plot(Time,Waves,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2);
xlabel('Time (s)'); ylabel('Units (A or V)'); title('Channel 1');

figure; plot(Time,Command,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2);
xlabel('Time (s)'); ylabel('Units (A or V)'); title('Channel 2');