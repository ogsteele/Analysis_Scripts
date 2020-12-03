%% NMDAR Subtraction Anlaysis
% Author: O.G. Steele
% Date: 09.10.20
% Description: 
%       Simple MATLAB script to perform subtraction analysis of the NMDAR
%       mediated component from compound events following the addition of
%       a selective NMDAR antagonist

%% Load in Data
% Use ephysIO.mat to load in the extracted ensembles from Eventer

% Select compound events first
title_str = "1. Select the ensemble average trace for compound events";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[file,path,~] = uigetfile('*.phy');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    % Navigate to directory and load file
    cd(path)
    COMP = ephysIO(file);
    clear('file','path','indx')

    % Select AMPAR events next
    title_str = "2. Select the ensemble average trace for AMPAR events";
    if ~ispc; menu(title_str,'OK'); end
    clear('title_str')
    [file,path,~] = uigetfile('*.phy');
    % Display file selection selection
    if isequal(file,0)
        disp('User selected Cancel')
        % If user selects cancel here, script will end here.
         return
    else
    disp(['User selected ', fullfile(path, file)])
    % Navigate to directory and load file
    cd(path)
    AMPA = ephysIO(file);
    clear('file','path','indx')
    end
end

%% Overlay the plots
% Overlay the plots of compound, AMPAR and calculated NMDAR

% Zero the ensembles by average the first 100 points
    % Compound
    basemean_COMP = mean(COMP.array(1:100,2));
    adjusted_COMP = COMP.array(:,2) - basemean_COMP;
    clear('basemean_COMP')
    % AMPAR
    basemean_AMPA = mean(AMPA.array(1:100,2));
    adjusted_AMPA = AMPA.array(:,2) - basemean_AMPA;
    clear('basemean_AMPA') 
    % Compensate for access resistance change
    %if exist('shift_percent','var') == 1
    %    adjusted_AMPA = adjusted_AMPA*((100 + shift_percent)/100);
    %else
    %end
% Calculate the NMDAR mediated trace
NMDA = AMPA;
NMDA.array(:,2) = adjusted_COMP - adjusted_AMPA;

% overlay the traces
    subfig = figure;
    % generate time as the size of the default event window in Eventer
    time = linspace(-0.01,0.04,length(adjusted_COMP));
    % convert amplitude from A to pA when plotting, raw still stored in Amps
    plot(time,(adjusted_COMP*(10^12)),'Linewidth',2)
    hold on
    plot(time,(adjusted_AMPA*(10^12)),'Linewidth',2)
    plot(time,(NMDA.array(:,2)*(10^12)),'Linewidth',2)
    % label graph
    xlabel('Time (Seconds)')
    ylabel('Amplitude (pA)')
    legend('Compound Event','+ 10 uM L-689,560 (AMPAR)','Calculated NMDAR current','Location','southeast','Linewidth',1)
    % stylise graph
    box off
    set(gca,'linewidth',2)
    set(gcf,'color','w');

% Tidy up
clear('ans','time','adjusted_AMPA','adjusted_COMP')
%  Navigate out of eventer's subfolders, to the root of the cell
cd ../../..
mkdir Subtraction_Results
save('Subtraction_Results/ensembles.mat','AMPA','COMP','NMDA');
saveas(subfig,'Subtraction_Results/sub_fig.pdf')