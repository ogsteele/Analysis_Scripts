%% Baseline Extractor


% Navigate the group directory and extract list of baseline.txts
path = cd(uigetdir()); % cd to ui choice (best to do by group at this point)
disp(['User selected ', path])
disp(['Parsing baseline.txt files ...'])
d = dir('**/*_baseline.txt'); % dir list of baseline files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden


% extract and concatenate into arrays
base = NaN(250,size(d,1));
for i = 1:size(d,1)
    a = table2cell(readtable(fullfile(d(i).folder,d(i).name)));
    b = cell2mat(a(:,2));
    % overwrite the selected bits
    base(1:size(b,1),i) = b;
    % disp((i/size(d,1))*100) show percentage progress
end

% Save as baseline.mat
prompt = {'End of compound ... (NMDAR blockade)'};
def = {'Gender_Geno'};
dlgtitle = 'Save baseline data as ...';
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,def);
save([char(answer),'_baseline.mat'],'base')
disp(['Saved as ', char(answer),'_baseline.mat in working directory'])
plot(base)
% tidy up
clear
