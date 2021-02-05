%% Event_splitter
% splits events and saves them as .phy files for stimfit 

%% clear workspace

close all
clear

%% Select .mat to load
title_str = "1. Load the ml_out variable of interest";
menu(title_str,'OK');
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Load the ml_out variable of interest');
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

%% calculate the NMDAR recordings
% subtract AMPAR from compound
NMDAR = ml_out.Compound.median - ml_out.AMPAR.median;
plot(NMDAR)
AMPAR = ml_out.AMPAR.median;
hold on
plot(AMPAR)
legend('NMDAR','AMPAR')
  
%% save recordings
% save as seperate .phy (use the filename variable)
time = 5e-5*[1:1001]';
% AMPAR
a = split(filename,'.');
name = append(char(a(1)),'_NMDAR.phy');
ephysIO(name,[time NMDAR],'s','A')
% NMDAR
a = split(filename,'.');
name = append(char(a(1)),'_AMPAR.phy');
ephysIO(name,[time AMPAR],'s','A')

% note: will need to concatenate at the end and then add in a time collumn
% as col 1
% x=5e-5*[1:1001]'
%% Notes for Oli as he's lazy
% process all of the waves and then concatenate into final structure
% a = load('ml_out.mat'); % will load a .mat as a the variable a
% plot(median([ml_out.Compound.median],2)); % will plot the median of the fields by row
% add in include or not row at the end ...
%   include = num2cell([1;0])
%    [test(:).data] = deal(include{:});