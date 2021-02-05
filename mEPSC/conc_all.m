%% Concatenate_all
% simple script to concatenate all contents of directory into one

%% Clear workspace

close all
clear
%% Select directory of .phy and list
title_str = "1. Navigate to directory for concatenation";
menu(title_str,'OK');
clear('title_str')
path = uigetdir();
% Display file selection selection
if isequal(path,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', path])
    % Navigate to directory and load file 
    cd(path)
end

% generate directory listing (and remove '.' files)
d=dir;
d=d(~ismember({d.name},{'.','..'}));
d = d(arrayfun(@(x) ~strcmp(x.name(1),'.'),d));

%% Load and concatenate

% Loop the loading and then concatenate
for i = 1:size(d,1)
    S = ephysIO(fullfile(d(i).folder, d(i).name));
    out(:,i) = S.array(:,2);
end

% add time as the first and save as ensemble average
time = 5e-5*[1:size(out,1)]';

% prompt user to name the output
prompt = {'Name the concatenated output'};
def = {'Genotype Gender Receptor'};
dlgtitle = 'output name';
dims = [1 50];
name = inputdlg(prompt,dlgtitle,dims,def);
% create the output in the structure
filename = append(name,'_ensemble_averages.phy');
% save the ensemble average
ephysIO(char(filename),[time out],'s','A') 