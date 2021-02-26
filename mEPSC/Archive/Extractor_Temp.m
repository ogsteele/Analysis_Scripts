% ExtractorTemp

% Ensure in cell directory

% Cell Information
    % split the path and copy to clipboard
    dirpath = split(pwd,'/');
    relevant = dirpath(12:17);
    % get the notes
    listing = dir('*notes.txt');
    notes = notesimport(fullfile(listing.folder,listing.name));
    relevant = [relevant; notes];
    mat2clip(relevant)
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the meta data yet?";
    menu(title_str,'Yup - proceed!');
    
    % whole cell capacitance    
    % select and load the wcp.mat file  
    listing = dir('*wcp.mat');                     
    load(fullfile(listing.folder,listing.name));   
    % extract the Cm value
    Cm = comp_output.comp_wcp.Cm;
    Cm = num2clip(Cm);
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the membrane capacitance data?";
    menu(title_str,'Well done - ');

% input resistance    
    % select and load the wcp.mat file  
    % listing = dir('*wcp.mat');                     % loaded above
    % load(fullfile(listing.folder,listing.name));   % loaded above
    % extract the Cm value
    Rin = comp_output.comp_wcp.Rin;
    Rin = num2clip(Rin);
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the input resistancee data?";
    menu(title_str,'Well done - ');
 
% membrane resistance    
    % select and load the wcp.mat file  
    % listing = dir('*wcp.mat');                     % loaded above
    % load(fullfile(listing.folder,listing.name));   % loaded above
    % extract the Cm value
    Rm = comp_output.comp_wcp.Rm;
    Rm = num2clip(Rm);
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the membrane resistancee data?";
    menu(title_str,'Well done - '); 
    
clear

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