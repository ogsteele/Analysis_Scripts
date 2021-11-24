function [Rheobase_pA,IR_MOhm,n_APs,dir_contents] = Rheobase_Wrapper

pA = [-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]';

pA_text = [{'-10 pA'},{'0 pA'},{'+10 pA'},{'+20 pA'},{'+30 pA'},{'+40 pA'},{'+50 pA'},...
    {'+60 pA'},{'+70 pA'},{'+80 pA'},{'+90 pA'},{'+100 pA'},{'+110 pA'},{'+120 pA'},...
    {'+130 pA'},{'+140 pA'},{'+150 pA'},{'+160 pA'},{'+170 pA'},{'+180 pA'}]';

% Get directory and list contents, removing empties (first two entries)
listing = dir(uigetdir);
dir_contents = listing(4:end,:);

% Preallocations
Rheobase_pA = zeros(1,size(dir_contents,1));
IR_MOhm = zeros(1,size(dir_contents,1));
n_APs = zeros(20,size(dir_contents,1));

    for i = 1:size(dir_contents,1)
        filename = fullfile(dir_contents(i).folder,dir_contents(i).name);
        [Rheobase_pA(i),IR_MOhm(i),n_APs(:,i)] = CurrentStep_Analysis(filename);
    end
    
end