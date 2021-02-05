%% Plotter
% simple script to plot overlays together in a nice presentation style

% clean environment and close figures
clear 
close all

% prompt user to name the output
prompt = {'Name the output'};
def = {'Genotype Gender'};
dlgtitle = 'output name';
dims = [1 50];
name = inputdlg(prompt,dlgtitle,dims,def);

%% Select AMPAR .h5 to load and plot
title_str = "1. load the AMPAR trace";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Load the AMPAR');
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
    AMPAR = ephysIO(file);
end

% plot the AMPAR trace
time = 5e-5*[1:1001]';
figure; plot(time, AMPAR.array(:,2),'linewidth',1)
set(gca,'linewidth',3,'fontsize',14)
set(gcf,'color','w');
box off
title(name)
ylabel('Amplitude (A)')
xlabel('Time (s)')

%% Select NMDAR .h5 to load
title_str = "1. load the NMDAR trace";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Load the NMDAR');
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
    NMDAR = ephysIO(file);
end

% overlay the NMDAR
hold on
plot(time, NMDAR.array(:,2),'linewidth',1)
legend1 = legend('AMPAR','NMDAR');
set(legend1,'LineWidth',1);

% save the output as pdf
cd ..
plotname = append(char(name),'_average_overlay.pdf');
saveas(gcf,plotname)