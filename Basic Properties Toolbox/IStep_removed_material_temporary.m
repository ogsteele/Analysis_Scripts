%% WrapperFunction

% check to see if I_step has already been run
try
    list = dir('**\*.mat'); % ending in .mat
    pat = 'CC_IV_1_full_';
    for i = 1:size(list,1)
        list(i).isfile = startsWith(list(i).name,pat); % is it one of ours?
    end
    fileFlags = [list.isfile];
    list = list(fileFlags); % save only the ones we're interested in
    istep_run = 1;
    warning('IStep run detected, running with extracted settings as defaults');
    load(list.name)
catch
    warning('Assuming IStep has not been run, running IStep normally');
    istep_run = 0;
end

% bridge balance corrected 
end
    elseif istep_run == 1
    balanced = output.Online_BB_performed;
    Rs_Init = output.Rs_Init;
    Vm_adjust = output.Offline_BB;
    IR = output.IR;
    Threshold = output.thresh;
    Afterhyperpolarisation = output.afterhyp;
    Overshoot = output.peak;
    offline_BB_performed = output.Offline_BB_performed;
    end



    
% What Genotype was the animal?
dlgTitle    = 'Genotpye';
dlgQuestion = 'What Genotype was this animal?';
if istep_run == 0
    genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4','APOE3');
elseif istep_run == 1
    genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4',(output.genotype));  
end
% Was the recording exposed to LPS/PBS?
dlgTitle    = 'Condition';
dlgQuestion = 'Which condition was used to generate this data?';
if istep_run == 0
    condition = questdlg(dlgQuestion,dlgTitle,'LPS (Test)','PBS (Control)','PBS (Control');
elseif istep_run == 1
    condition = questdlg(dlgQuestion,dlgTitle,'LPS (Test)','PBS (Control)',(output.condition));
end

% What was the slice ID, and  add any notes?
prompt = {'What was the slice ID for this recording?','Would you like to add any notes?'};
dlgtitle = 'Slice ID & notes';
if istep_run == 0
    definput = {'eg. E420211203#1','eg. dodgy input'};
elseif istep_run == 1
    definput = {char(output.ID),char(output.Notes)};
end
dims = [1 40];
slice_id_notes = inputdlg(prompt,dlgtitle,dims,definput);
