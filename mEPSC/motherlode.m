%% Tidy up 

clear
close all

%% Navigation 
% Navigate to the folder containing all analysed data
path = cd(uigetdir()); % cd to ui choice (best to do by group at this point)
disp(['User selected ', path])

%% ml_out 

% begin parsing ml_out.mat files for event information
disp(['Parsing ml_out.mat files ...'])
tic
d = dir('**/ml_out.mat'); % dir list of wcp.mat files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden
toc

% extract the desired information and parse into a structure
numIterations  = size(d,1);
ppm = ParforProgressbar(numIterations,...
    'title','ml out.mat extraction');
parfor i = 1:size(d,1)
    a = load(fullfile(d(i).folder,d(i).name));
    % create ID
    split_folder = split(d(i).folder,'/');
    out(i).Genotype = char(split_folder(end-5));
    out(i).Sex = char(split_folder(end-4));
    out(i).ID = split_folder(end-2);
    out(i).Slice = split_folder(end-1);
    out(i).Cell = split_folder(end);
    % extract AMPAR values
    out(i).AMPA_event_num = a.ml_out.AMPAR.event_num;
    out(i).AMPA_Hz = a.ml_out.AMPAR.event_num/600;
    out(i).AMPA_Amp = a.ml_out.AMPAR.event_amp;
    % extract Compound values
    out(i).Comp_event_num = a.ml_out.Compound.event_num;
    out(i).Comp_Hz = a.ml_out.Compound.event_num/600;
    out(i).Comp_Amp = a.ml_out.Compound.event_amp;
    ppm.increment();   
end
delete(ppm); clear ppm;

%% wcp.mat 

disp(['Parsing wcp.mat files ...'])
tic
d = dir('**/*_wcp.mat'); % dir list of wcp.mat files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden
toc

% extract the desired information and parse into a structure
raw_Rs = zeros(180,size(d,1));
comp_Rs = zeros(180,size(d,1));

numIterations  = size(d,1);
ppm = ParforProgressbar(numIterations,...
    'title','wcp out.mat extraction');
parfor i = 1:size(d,1)
    a = load(fullfile(d(i).folder,d(i).name));
    out(i).raw_Rs = a.comp_output.raw_wcp.Rs(1:180);
    out(i).mean_raw_Rs = mean(a.comp_output.raw_wcp.Rs(1:180));
    out(i).comp_Rs = a.comp_output.comp_wcp.Rs(1:180);
    out(i).mean_comp_Rs = mean(a.comp_output.comp_wcp.Rs(1:180));
    ppm.increment();   
end
delete(ppm); clear ppm;

%% Splitting for later convenience

E4_subset = out(strcmp({out.Genotype}, 'APOE4'));
E3_subset = out(strcmp({out.Genotype}, 'APOE3'));

M_subset = out(strcmp({out.Sex}, 'Male'));
F_subset = out(strcmp({out.Sex}, 'Female'));

E4M_subset = out(strcmp({out.Genotype}, 'APOE4') & strcmp({out.Sex}, 'Male'));
E4F_subset = out(strcmp({out.Genotype}, 'APOE4') & strcmp({out.Sex}, 'Female'));

E3M_subset = out(strcmp({out.Genotype}, 'APOE3') & strcmp({out.Sex}, 'Male'));
E3F_subset = out(strcmp({out.Genotype}, 'APOE3') & strcmp({out.Sex}, 'Female'));
%% Plotting

% plot Rs values across all conditions
x = linspace(0,1,180);
figure; shadedErrorBar(x,mean([out.raw_Rs],2),std([out.raw_Rs],[],2),'lineprops','r'); 
hold on; shadedErrorBar(x,mean([out.comp_Rs],2),std([out.comp_Rs],[],2),'lineprops','b');
xlabel('sweep'); ylabel('Rs (MOhm)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('raw','comp')

% plot event frequency between condition
