%% Baseline Extractor


% Navigate the group directory and extract list of baseline.txts
path = cd(uigetdir()); % cd to ui choice (best to do by group at this point)
disp(['User selected ', path])
disp(['Parsing *_rscomp.phy files ...'])
 disp(['Progress ...'])
    disp('0.0000') %show percentage progress
d = dir('**/*_rscomp.phy'); % dir list of baseline files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden

% extract and concatenate into arrays
base = NaN(250,size(d,1)); % recordings shouldn't be longer than 250 waves
for ii = 1:size(d,1)
    S = ephysIO(fullfile(d(ii).folder,d(ii).name));
    % create a matrix of the split recording and preallocate
    % create list of start and end points determined by the sample rate, length
    % of split, length of recording
    n_splits = floor(((size(S.array(:,2),1))/200000));
    length_raw_split = 200000;
    % round down n_splits to not deal with incomplete sections and preallocate
    start(:,1) = zeros(n_splits,1);
    finish(:,1) = zeros(n_splits,1);
    
    % define the starting numbers
    start(1,1) = 1;
    finish(1,1) = 1 + length_raw_split-1;
    for i = 2:n_splits
        start(i,1) = finish(i-1,1) + 1;
        finish(i,1) = start(i,1) + length_raw_split-1;
    end

    % create a matrix of the split recording and preallocate
    splits = zeros(length_raw_split,n_splits);
    for s = 1:n_splits
        splits(:,s) = S.array(start(s):finish(s),2);
    end
    
    clear('start','finish')
    % median of each split
    b = median(splits,1)';
    % overwrite the selected bits
    base(1:size(b,1),ii) = b;
    disp(['Progress ...'])
    disp((ii/size(d,1))*100) %show percentage progress
end
disp(['Complete'])

% Save as baseline.mat
prompt = {'End of compound ... (NMDAR blockade)'};
def = {'Gender_Geno'};
dlgtitle = 'Save baseline data as ...';
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,def);
figure
plot(base)
save([char(answer),'_baseline.mat'],'base')
disp(['Saved as ', char(answer),'_baseline.mat in working directory'])

% tidy up
clear
