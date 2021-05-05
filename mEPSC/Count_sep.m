%% Count_sep

% script to count the events, and then seperate dependant on region.
% this provides you with Compund and AMPAR events & medians

% will also save compound, AMPAR and NMDAR seperated sections as `ml_out`

% clear workspace to begin
close all
clear

%% Load in w/ readtable
% Select event counts to visualise
title_str = "1. Select eventer analysis directory (above eventer.output)";
menu(title_str,'OK');
clear('title_str')

filepath = uigetdir();
%[file,path,~] = uigetfile('*.*','1. Select event counts');
% Display file selection selectio
%if isequal(file,0)
if isequal(filepath,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', filepath])
    %disp(['User selected ', fullfile(path, file)])
    filename = filepath;
    %filename = file;
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

disp(['Total number of waves = ', num2str(waves)])

%% Wave locations
% ask for the AMPAR / NMDAR / Compound wave times 
% necessary in case a recording is short
prompt = {'End of compound ... (NMDAR blockade)',...
    'Start of AMPAR ... (Set to 0 if NA)'};
def = {'60','120'};
dlgtitle = 'Wave numbers';
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,def);
answer = str2double(answer); % convert to number
Comp_end = answer(1);
AMPAR_start = answer(2);

% number of events at each point
cum_counts = cumsum(counts);
comp_event_num = cum_counts(Comp_end);

%% Load in w/ ePhysIO
% Select event data to visualise
% already in correct directory
S = ephysIO(fullfile(filepath,'eventer.output/All_events/event_data.phy'));

%% Seperate and average sections
% according to drug application and wave

Comp_events = S.array(:,2:comp_event_num);
Comp_median = median(Comp_events,2);

if AMPAR_start < Comp_end
    disp('No AMPAR isolated events detected')
    figure
    plot(Comp_median)
    legend('Compound')
    title('Event overlays')
    ylabel('Amplitude (A)')
    xlabel('Data points')
else
    AMPAR_event_num = cum_counts(AMPAR_start);
    AMPAR_events = S.array(:,AMPAR_event_num:end);
    AMPAR_median = median(AMPAR_events,2);
    figure
    plot(Comp_median)
    hold on
    plot(AMPAR_median)
    plot(Comp_median - AMPAR_median)
    legend('Compound','AMPAR','NMDAR')
    title('Event overlays')
    ylabel('Amplitude (pA)')
    xlabel('Data points')
end

%% save outputs
% navigate out of ml_out directory so data is not overridden
cd ../../..

% Cell Information
    % split the path and copy to clipboard
    dirpath = split(pwd,'/');
    relevent = dirpath(end-5:end);
    % get the notes
    listing = dir('*notes.txt');
    listing = listing(~startsWith({listing.name}, '.')); % remove deleted/hidden
    notes = notesimport(fullfile(listing.folder,listing.name));
    
ml_out.Info.Genotype = relevent(1);
ml_out.Info.Gender = relevent(2);
ml_out.Info.Date = relevent(3);
ml_out.Info.Animal_Num = relevent(4);
ml_out.Info.Slice_Num = relevent(5);
ml_out.Info.Cell_Num = relevent(6);
ml_out.Info.Notes = notes;

ml_out.All_Raw.events = S.array(:,2:end);

ml_out.Compound.events = Comp_events;
ml_out.Compound.median = Comp_median;
ml_out.Compound.event_num = size(Comp_events,2);
ml_out.Compound.event_amp = min(Comp_median)*1e12;

if AMPAR_start > Comp_end
    ml_out.AMPAR.events = AMPAR_events;
    ml_out.AMPAR.median = AMPAR_median;
    ml_out.AMPAR.event_num = size(AMPAR_events,2);
    ml_out.AMPAR.event_amp = min(AMPAR_median)*1e12;
else
end

save('ml_out.mat','ml_out')


%% Do the job of legacy script `event_sep`
NMDAR = ml_out.Compound.median - ml_out.AMPAR.median;
AMPAR = ml_out.AMPAR.median;
COMP = ml_out.Compound.median;
% save recordings
% save as seperate .phy (use the filename variable)
time = 5e-5*[1:size(AMPAR,1)]';

% NMDAR
a = split(filename,'/');
name = append(char(a(end-3)),'_',char(a(end-2)),'_',char(a(end-1)),'_average_NMDAR.phy');
%name = append(char(a(1)),'_NMDAR.phy');
ephysIO(name,[time NMDAR],'s','A')

% AMPAR
name = append(char(a(end-3)),'_',char(a(end-2)),'_',char(a(end-1)),'_average_AMPAR.phy');
%name = append(char(a(1)),'_AMPAR.phy');
ephysIO(name,[time AMPAR],'s','A')

% Compound
name = append(char(a(end-3)),'_',char(a(end-2)),'_',char(a(end-1)),'_average_COMP.phy');
%name = append(char(a(1)),'_AMPAR.phy');
ephysIO(name,[time COMP],'s','A')
%%

function [notes1] = notesimport(filename) 
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, 4];
opts.Delimiter = "";

% Specify column names and types
opts.VariableNames = "UsableForCompoundEventAnalysis";
opts.VariableTypes = "char";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "UsableForCompoundEventAnalysis", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "UsableForCompoundEventAnalysis", "EmptyFieldRule", "auto");

% Import the data
notes1 = readtable(filename, opts);

% Convert to output type
notes1 = table2cell(notes1);
numIdx = cellfun(@(x) ~isnan(str2double(x)), notes1);
notes1(numIdx) = cellfun(@(x) {str2double(x)}, notes1(numIdx));

% Clear temporary variables
clear opts
end