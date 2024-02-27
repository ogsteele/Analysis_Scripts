[file, path] = uigetfile('*.*');
S = ephysIO(fullfile(path,file));
figure; plot(S.array(:,1),S.array(:,2:end))
% ask is this exemplary
dlgTitle    = 'Exemplary?';
dlgQuestion = 'Is this recording considered exemplary?';
run = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
if run == "Yes"
    % if run = yes, go ahead and ask the user to pick an exemplary wave
    disp('enter current value (pA) to analyse spike frequency adaptation at')
    prompt = {...
        'Enter adaptation stimulation value (pA):'};
    dlg_title = 'Adaptation stimulation value input';
    num_lines = 1;
    def = {'300'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
    answer = str2double(answer);
    pA = [-200:20:400];
    step = find(pA==300);
    figure; plot(S.array(:,1), S.array(:,step)*1000)
    xlabel('Time (s)')
    ylabel('Membrane Potential (mV)')
    % would you like to save this as an emf? 
    dlgTitle    = 'Save as emf?';
    dlgQuestion = 'Do you want a figure out of this trace?';
    emf = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
    if emf == "Yes"
        saveas(gcf,'exemplary.emf')
    else 
    end
else
end