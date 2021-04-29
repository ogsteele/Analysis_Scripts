%% wcp_extractor

% Navigate the group directory and extract list of wcp.mat files
path = cd(uigetdir()); % cd to ui choice (best to do by group at this point)
disp(['User selected ', path])
disp(['Parsing wcp.mat files ...'])
d = dir('**/*_wcp.mat'); % dir list of wcp.mat files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden

% extract the desired information and parse into a structure
raw_Rs = zeros(180,size(d,1));
comp_Rs = zeros(180,size(d,1));
parfor i = 1:size(d,1)
    a = load(fullfile(d(i).folder,d(i).name));
    raw_Rs(:,i) = a.comp_output.raw_wcp.Rs(1:180);
    comp_Rs(:,i) = a.comp_output.comp_wcp.Rs(1:180);
end
    
% plot Rs values
x = linspace(0,1,180);
figure; shadedErrorBar(x,mean(raw_Rs,2),std(raw_Rs,[],2),'lineprops','r'); 
hold on; shadedErrorBar(x,mean(comp_Rs,2),std(comp_Rs,[],2),'lineprops','b');
xlabel('sweep'); ylabel('Rs (MOhm)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('raw','comp')
    
    





%% NOTES
% assign variables to an output structure
out.Ih = Ih; % holding current in pA
out.Ip = Ip; % peak amplitude of the negative current transient in pA
out.base = base; % pre-pulse baseline in pA
out.peak = peak; % peak
out.Rs = Rs; % series resistance in Mohm
out.Iss = Iss; % intrapulse steady state current in pA
out.Rin = Rin; % input resistance in Mohm
out.Rm = Rm; % specific membrane resistance in Mohm
out.Cm = Cm; % capacitance in pF
out.Q = Q; % charge in C