%% mEPSP_Sub
%
% Author: O.G. Steele
% Date (initial): 29.05.21

% Description:
%   Simple code to extract averages of a subtraction experiment from
%   Eventer.output based on the time of drug application. 

% Updates: 
%
%
%
%% Get amplifier parameters
% Read in amplifier settings from notes.txt

% split length in seconds
split_length = 10;

%% Time points
% ask the user the time points of the experiment
prompt = {...
    'Start of Recording (s)',...
    'Application of L-689,560 (s)',...
    'End of Recording (s)',...
    };
defaults = {'0','600','1800'};
dlgtitle = 'Time points in recording';
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,defaults);
answer = str2double(answer); % convert to number

% if 'answer' is empty, terminate script
if isempty(answer)
    disp('User cancelled on time point entry')
    return;
else
    Start_s = answer(1); % seconds
    Start_w = Start_s / split_length; % wave
    Drug_s = answer(2); % seconds
    Drug_w = Drug_s / split_length; % wave
    End_s = answer(3); % seconds
    End_w = End_s / split_length; % wave
end

%% Plot cumulative sum of events
% Select event counts to visualise
title_str = "1. Select eventer analysis directory (above eventer.output)";
menu(title_str,'OK');
clear('title_str')

filepath = uigetdir();
if isequal(filepath,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', filepath])
    filename = filepath;
    % Navigate to directory and load file with ephysIO
    count_filepath = fullfile(filepath,'eventer.output/All_events/event_counts.txt');
    cd(fullfile(filepath,'eventer.output/All_events'));
    counts = table2array(readtable(count_filepath));
    waves = size(counts,1);
end
% tidy workspace
figure
clear('path','ans')
plot(cumsum(counts))
title('Cumulative sum of events')
xlabel('Wave number')
ylabel('Event number')

% display summary counts
disp(['Total number of waves = ', num2str(waves)])
disp(['Total number of events = ', num2str(sum(counts))])
disp(['Total average of frequency = ', num2str(sum(counts)/(waves*split_length))])


%% Load in w/ ePhysIO
% Select event data to visualise
% already in correct directory
S = ephysIO(fullfile(filepath,'eventer.output/All_events/event_data.phy'));

%% Seperate and average sections
% according to drug application and wave
% taken from Count_Sep.m, but needs to be made applicable to this set up. 

% Comp_events = S.array(:,2:comp_event_num);
% Comp_median = median(Comp_events,2);

% Define points of interest
% S.array(:,1) = time
Compound_events = S.array(:,Start_w+2:Drug_w+2);
AMPAR_events = S.array(:,Drug_w+2:End_w);

comp = median(
ampa

baseline = table2array(readtable(d.name));
t = comp.array(:,1);
b = median(baseline(10:80,2)); % 'b' is the average baseline before L689560
a = median(baseline(end-70:end,2)); % 'a' is the average baseline after L689560
A = (a + ampa.array(:,2)) * b/a; % corrected AMPAR
B = b + comp.array(:,2); % corrected compound
C = (b + (B-A)); % corrected NMNDAR

figure
hold on
plot(median(Compound_events,2))
plot(median(AMPAR_events,2))
legend('Comp','AMPAR')
plot((median(Compound_events,2))-(median(AMPAR_events,2)))
legend('Baseline','+ L-689,560','Difference')
title('Event overlays')
ylabel('Amplitude (V)')
xlabel('Data points')
box off
set(gcf,'color','white')
set(gca,'linewidth',2)
hold off


% if AMPAR_start < Comp_end
%     disp('No AMPAR isolated events detected')
%     figure
%     plot(Comp_median)
%     legend('Compound')
%     title('Event overlays')
%     ylabel('Amplitude (A)')
%     xlabel('Data points')
% else
%     AMPAR_event_num = cum_counts(AMPAR_start);
%     AMPAR_events = S.array(:,AMPAR_event_num:end);
%     AMPAR_median = median(AMPAR_events,2);
%     figure
%     plot(Comp_median)
%     hold on
%     plot(AMPAR_median)
%     plot(Comp_median - AMPAR_median)
%     legend('Compound','AMPAR','NMDAR')
%     title('Event overlays')
%     ylabel('Amplitude (pA)')
%     xlabel('Data points')
% end