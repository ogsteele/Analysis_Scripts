%% overlay_plotter
% simple script that can plot the median events and then in a
% slightly darker shade plot the median overlay.


% clean environment and close figures
clear 
close all

% prompt user to name the output
prompt = {'Name the output'};
def = {'Genotype Gender'};
dlgtitle = 'output name';
dims = [1 50];
name = inputdlg(prompt,dlgtitle,dims,def);

%% Select all AMPAR .phy to load and plot
title_str = "1. Select all AMPAR .phy to load and plot";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','Select all AMPAR .phy to load and plot');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    filename = file;
    % Navigate to directory and load file 
    cd(path)
    all = ephysIO(file);
end

% plot the AMPAR trace
time = 5e-5*[1:1001]';
%figure; plot(time, all.array(:,2:end),'linewidth',0.5,'color',[0.8500, 0.3250, 0.0980, 0.3]) % orange
% figure; plot(time, all.array(:,2:end),'linewidth',0.5,'color',[0, 0.4470,0.7410, 0.3]) % blue
 figure; plot(time, all.array(:,2:end),'linewidth',0.5,'color',[0.9290,0.6940, 0.1250 0.3]) % yellow
set(gca,'linewidth',3,'fontsize',14)
set(gcf,'color','w');
box off
ylabel('Amplitude (A)')
xlabel('Time (s)')
ylim([-1e-11 0.5e-11])

%% Select average AMPAR .h5 to load and plot
title_str = "Select average AMPAR .h5 to load and plot";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','Select average AMPAR .h5 to load and plot');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    filename = file;
    % Navigate to directory and load file 
    cd(path)
    average = ephysIO(file);
end

hold on
% plot(time, average.array(:,2),'linewidth',2,'color',[0.8500, 0.3250, 0.0980, 1]) % orange
% plot(time, average.array(:,2),'linewidth',2,'color',[0, 0.4470, 0.7410, 1]) % blue
 plot(time, average.array(:,2),'linewidth',2,'color',[0.9290, 0.6940, 0.1250, 1]) % yellow
set(gca,'visible','off')

% save fig
plotname = append(char(name),'_events.pdf');
saveas(gcf,plotname)
%%
%S = ephysIO('average_average.h5')
%s = ephysIO('E3F_AMPAR_ensemble_averages.phy')
%s = ephysIO('E3F_AMPAR_ensemble_averages.phy')
%figure
%plot(s.array(:,2:end),'color',[1,0,0,0.2])
%plot(S.array(:,2),'color',[1,0,0,1],'linewidth',2)
%hold on
%plot(s.array(:,2:end),'color',[1,0,0,0.2])