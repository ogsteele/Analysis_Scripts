% Ensure in the cell directory
% Simple script to help extract data from the storage created with the the
% other mEPSC scripts around here

% Cell Information
    % split the path and copy to clipboard
    dirpath = split(pwd,'/');
    relevant = dirpath(10:15);
    % get the notes
    listing = dir('*notes.txt');
    notes = notesimport(fullfile(listing.folder,listing.name));
    relevant = [relevant; notes];
    mat2clip(relevant)
    
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the meta data yet?";
    menu(title_str,'Yup - proceed!');
        
% Noise levels
    % select *noise.mat file and load
    listing = dir('*noise.mat');
    load(fullfile(listing.folder,listing.name));
    % convert noise levels to array
    noise_array = table2array(noise_output.noise_levels);
    % copy to clipboard
    noise_array = num2clip(noise_array(:,2));
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the noise data?";
    menu(title_str,'Yep, next set');
        
% Time constants    
    % select and load the summary.txt file
    listing = dir('*eventer.output/ALL_events/summary.txt');
    % organise table for easier extraction
    sum = readtable(fullfile(listing.folder,listing.name));
    sum = table(sum.Var2,'RowNames',sum.Var1);
    % pull out values
    rise = table2array(sum({'Rise time constant of the model PSC fit (ms)'},:));
    decay = table2array(sum({'Decay time constant of the model PSC fit (ms)'},:));
    rise_decay = [rise;decay];
    rise_decay = num2clip(rise_decay);
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the time constant data?";
    menu(title_str,'Duhh ... ');
    
% Tonic current    
    % select and load the wcp.mat file  
    listing = dir('*wcp.mat');
    load(fullfile(listing.folder,listing.name));
    % extract the Ih value
    Ih = comp_output.raw_wcp.Ih;
    Ih = num2clip(Ih);
    % PASTE NOW INTO EXCEL
    title_str = "Have you pasted the Tonic Current data?";
    menu(title_str,'Well done - ');
    
% whole cell capacitance    
    % select and load the wcp.mat file  
    % listing = dir('*wcp.mat');                     % loaded above
    % load(fullfile(listing.folder,listing.name));   % loaded above
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

function out = mat2clip(a, delim)

%MAT2CLIP  Copies matrix to system clipboard.
%
% MAT2CLIP(A) copies the contents of 2-D matrix A to the system clipboard.
% A can be a numeric array (floats, integers, logicals), character array,
% or a cell array. The cell array can have mixture of data types.
%
% Each element of the matrix will be separated by tabs, and each row will
% be separated by a NEWLINE character. For numeric elements, it tries to
% preserve the current FORMAT. The copied matrix can be pasted into
% spreadsheets.
%
% OUT = MAT2CLIP(A) returns the actual string that was copied to the
% clipboard.
%
% MAT2CLIP(A, DELIM) uses DELIM as the delimiter between columns. The
% default is tab (\t).
%
% Example:
%   format long g
%   a = {'hello', 123;pi, 'bye'}
%   mat2clip(a);
%   % paste into a spreadsheet
%
%   format short
%   data = {
%     'YPL-320', 'Male',   38, true,  uint8(176);
%     'GLI-532', 'Male',   43, false, uint8(163);
%     'PNI-258', 'Female', 38, true,  uint8(131);
%     'MIJ-579', 'Female', 40, false, uint8(133) }
%   mat2clip(data);
%   % paste into a spreadsheet
%
%   mat2clip(data, '|');   % using | as delimiter
%
% See also CLIPBOARD.

% VERSIONS:
%   v1.0 - First version
%   v1.1 - Now works with all numeric data types. Added option to specify
%          delimiter character.
%
% Copyright 2009 The MathWorks, Inc.
%
% Inspired by NUM2CLIP by Grigor Browning (File ID: 8472) Matlab FEX.

error(nargchk(1, 2, nargin, 'struct'));

if ndims(a) ~= 2
  error('mat2clip:Only2D', 'Only 2-D matrices are allowed.');
end

% each element is separated by tabs and each row is separated by a NEWLINE
% character.
sep = {'\t', '\n', ''};

if nargin == 2
  if ischar(delim)
    sep{1} = delim;
  else
    error('mat2clip:CharacterDelimiter', ...
      'Only character array for delimiters');
  end
end

% try to determine the format of the numeric elements.
switch get(0, 'Format')
  case 'short'
    fmt = {'%s', '%0.5f' , '%d'};
  case 'shortE'
    fmt = {'%s', '%0.5e' , '%d'};
  case 'shortG'
    fmt = {'%s', '%0.5g' , '%d'};
  case 'long'
    fmt = {'%s', '%0.15f', '%d'};
  case 'longE'
    fmt = {'%s', '%0.15e', '%d'};
  case 'longG'
    fmt = {'%s', '%0.15g', '%d'};
  otherwise
    fmt = {'%s', '%0.5f' , '%d'};
end

if iscell(a)  % cell array
    a = a';
    
    floattypes = cellfun(@isfloat, a);
    inttypes = cellfun(@isinteger, a);
    logicaltypes = cellfun(@islogical, a);
    strtypes = cellfun(@ischar, a);
    
    classType = zeros(size(a));
    classType(strtypes) = 1;
    classType(floattypes) = 2;
    classType(inttypes) = 3;
    classType(logicaltypes) = 3;
    if any(~classType(:))
      error('mat2clip:InvalidDataTypeInCell', ...
        ['Invalid data type in the cell array. ', ...
        'Only strings and numeric data types are allowed.']);
    end
    sepType = ones(size(a));
    sepType(end, :) = 2; sepType(end) = 3;
    tmp = [fmt(classType(:));sep(sepType(:))];
    
    b=sprintf(sprintf('%s%s', tmp{:}), a{:});

elseif isfloat(a)  % floating point number
    a = a';
    
    classType = repmat(2, size(a));
    sepType = ones(size(a));
    sepType(end, :) = 2; sepType(end) = 3;
    tmp = [fmt(classType(:));sep(sepType(:))];
    
    b=sprintf(sprintf('%s%s', tmp{:}), a(:));

elseif isinteger(a) || islogical(a)  % integer types and logical
    a = a';
    
    classType = repmat(3, size(a));
    sepType = ones(size(a));
    sepType(end, :) = 2; sepType(end) = 3;
    tmp = [fmt(classType(:));sep(sepType(:))];
    
    b=sprintf(sprintf('%s%s', tmp{:}), a(:));

elseif ischar(a)  % character array
    % if multiple rows, convert to a single line with line breaks
    if size(a, 1) > 1
      b = cellstr(a);
      b = [sprintf('%s\n', b{1:end-1}), b{end}];
    else
      b = a;
    end
    
else
    error('mat2clip:InvalidDataType', ...
      ['Invalid data type. ', ...
      'Only cells, strings, and numeric data types are allowed.']);

end

clipboard('copy', b);

if nargout
  out = b;
end
end 

function arraystring = num2clip(array)
%NUM2CLIP copies a numerical-array to the clipboard
%   
%   ARRAYSTRING = NUM2CLIP(ARRAY)
%   
%   Copies the numerical array ARRAY to the clipboard as a tab-separated
%   string.  This format is suitable for direct pasting to Excel and other
%   programs.
%   
%   The tab-separated result is returned as ARRAYSTRING.  This
%   functionality has been included for completeness.
%   
%Author: Grigor Browning
%Last update: 02-Sept-2005

%convert the numerical array to a string array
%note that num2str pads the output array with space characters to account
%for differing numbers of digits in each index entry
arraystring = num2str(array); 
arraystring(:,end+1) = char(10); %add a carrige return to the end of each row
%reshape the array to a single line
%note that the reshape function reshape is column based so to reshape by
%rows one must use the inverse of the matrix
arraystring = reshape(arraystring',1,prod(size(arraystring))); %reshape the array to a single line

arraystringshift = [' ',arraystring]; %create a copy of arraystring shifted right by one space character
arraystring = [arraystring,' ']; %add a space to the end of arraystring to make it the same length as arraystringshift

%now remove the additional space charaters - keeping a single space
%charater after each 'numerical' entry
arraystring = arraystring((double(arraystring)~=32 | double(arraystringshift)~=32) & ~(double(arraystringshift==10) & double(arraystring)==32) );

arraystring(double(arraystring)==32) = char(9); %convert the space characters to tab characters

clipboard('copy',arraystring); %copy the result to the clipboard ready for pasting

end 

