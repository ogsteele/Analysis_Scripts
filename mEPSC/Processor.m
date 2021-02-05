%% Processing
% Simple script for processing and concatenating the dataset

%% Clear workspace
clear 
close all

%% Load in the ml_out variable
% Select event counts to visualise
title_str = "1. Load the ml_out variable of interest";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Select event counts');
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
    load(file);
end

%% Unique ID
% prompt user to enter the Unique ID as designated in note book
% GenoGender_Number_Rep, eg "E3F_897_1"
prompt = {'Unique ID'};
def = {'GenoGender_Number_Rep'};
dlgtitle = 'Unique ID';
dims = [1 50];
UID = inputdlg(prompt,dlgtitle,dims,def);
% create the output in the structure
ml_out.UID = string(UID);

%% Plot the data
% allows the user to see the data and make the notes accordingly
figure
plot(ml_out.Compound.median)
hold on
exist_ans = isfield(ml_out,'AMPAR');
    if exist_ans == 1
        plot(ml_out.AMPAR.median)
        plot(ml_out.Compound.median - ml_out.AMPAR.median)
        legend('Compound','AMPAR','NMDAR')
        title('Event overlays')
        ylabel('Amplitude (pA)')
        xlabel('Data points')
    else
        % create empty AMPAR only events if not included
        ml_out.AMPAR.events = [];
        ml_out.AMPAR.median = [];
        ml_out.AMPAR.event_num = [];
        ml_out.AMPAR.event_amp = [];
        % create original plot with only AMPAR
        legend('Compound')
        title('Event overlays')
        ylabel('Amplitude (pA)')
        xlabel('Data points')        
    end    

% save the plot to the current wd as a pdf
file = append(string(UID),'_subtractionPlot.pdf');
saveas(gcf,file)

%% User notes
% prompt user to enter helpful user notes for quick evaluation"
prompt = {'User notes'};
def = {'Compound Only etc ... '};
dlgtitle = 'User notes';
dims = [1 50];
usernotes = inputdlg(prompt,dlgtitle,dims,def);
% create the output in the structure
ml_out.user_notes = string(usernotes);

%% Save ml_out
% overwrite old ml_out.mat
name = append(string(UID),'_ml_out.mat');
save(name,'ml_out')

%% Notes for Oli as he's lazy
% process all of the waves and then concatenate into final structure
% a = load('ml_out.mat'); % will load a .mat as a the variable a
% plot(median([ml_out.Compound.median],2)); % will plot the median of the fields by row
% add in include or not row at the end ...
%   include = num2cell([1;0])
%    [test(:).data] = deal(include{:});