%% motherlode 

% requirements:
%    - parforprogress by Frerk Saxen (Mathworks File Exchange, 2019) & 
%    shaded_error by Rob Campbell (Mathworks File Exchange, 2018) 
%    (both available on mathworks, but also in the `mEPSC/third_party` 
%    github page)
%    - pargorprogress requires the instrument control toolbox
%    (https://www.mathworks.com/products/instrument.html) and possibly the
%    parralel processing toolbox
%
% not an interesting script, just a timesaving one. Will run through all
% the data and pull out interesting features as desired - so far ... 
%
%     create ID
%     out(i).Genotype = char(split_folder(end-5));
%     out(i).Sex = char(split_folder(end-4));
%     out(i).ID = split_folder(end-2);
%     out(i).Slice = split_folder(end-1);
%     out(i).Cell = split_folder(end);
%     extract AMPAR values
%     out(i).AMPA_event_num = a.ml_out.AMPAR.event_num;
%     out(i).AMPA_Hz = a.ml_out.AMPAR.event_num/600;
%     out(i).AMPA_Amp = a.ml_out.AMPAR.event_amp;
%     extract Compound values
%     out(i).Comp_event_num = a.ml_out.Compound.event_num;
%     out(i).Comp_Hz = a.ml_out.Compound.event_num/600;
%     out(i).Comp_Amp = a.ml_out.Compound.event_amp;
%     extract whole cell properties
%     out(i).raw_Rs = a.comp_output.raw_wcp.Rs(1:180);
%     out(i).mean_raw_Rs = mean(a.comp_output.raw_wcp.Rs(1:180));
%     out(i).comp_Rs = a.comp_output.comp_wcp.Rs(1:180);
%     out(i).mean_comp_Rs = mean(a.comp_output.comp_wcp.Rs(1:180));
%     extract event counts
%     out(i).event_counts = NaN(180,1);
%     out(i).cumsum_event_counts = cumsum(out(i).event_counts);
%     extract baseline data
%     out(i).baseline = a(1:180,2);
%     out(i).early_base = mean(out(i).baseline(1:60));
%     out(i).late_base = mean(out(i).baseline(121:180));
%     
%     motherlode will save the above into an output structure, whilst also subsetting the data and saving to the workspace only currently. 
%     motherlode will produce several plots of the data also.
%   
%     TODO:
%            include more wcp (ie, all of them) 
%            include more actual recording data
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


%% Cumsum Hz

disp(['Parsing event_counts.txt files ...'])
tic
d = dir('**/mlm_500_new/eventer.output/ALL_events/event_counts.txt'); % dir list of wcp.mat files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden
toc

numIterations  = size(d,1);
ppm = ParforProgressbar(numIterations,...
    'title','event_counts.txt extraction');
for i = 1:size(d,1)
    a = table2array(readtable(fullfile(d(i).folder,d(i).name)));
    out(i).event_counts = NaN(180,1);
    if size(a,1) >= 180
        out(i).event_counts = a(1:180);
    else 
        out(i).event_counts(1:size(a,1)) = a;
    end
    out(i).cumsum_event_counts = cumsum(out(i).event_counts);
    ppm.increment(); 
end
delete(ppm); clear ppm;


%% Baseline_txt

% Navigate the group directory and extract list of baseline.txts
disp(['Parsing baseline.txt files ...'])
tic
d = dir('**/*_baseline.txt'); % dir list of baseline files
d = d(~startsWith({d.name}, '.')); % remove deleted/hidden
toc

numIterations  = size(d,1);
ppm = ParforProgressbar(numIterations,...
    'title','baseline.txt extraction');
for i = 1:size(d,1)
    a = table2array(readtable(fullfile(d(i).folder,d(i).name)));
    out(i).baseline = NaN(180,1);
    if size(a,1) >= 180
        out(i).baseline = a(1:180,2);
    else 
        out(i).baseline(1:size(a,1)) = a(:,2);
    end
    out(i).early_base = mean(out(i).baseline(1:60));
    out(i).late_base = mean(out(i).baseline(121:180));
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

    %% plot Rs values across all conditions
x = linspace(0,180,180);
figure; shadedErrorBar(x,mean([out.raw_Rs],2),std([out.raw_Rs],[],2),'lineprops','r'); 
hold on; shadedErrorBar(x,mean([out.comp_Rs],2),std([out.comp_Rs],[],2),'lineprops','b');
xlabel('sweep'); ylabel('Rs (MOhm)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('raw','comp')

    %% plot event frequency between condition
    figure; 
X = categorical({'E3M','E3F','E4M','E4F'});
%X = reordercats(X,{'E3M','E3F','E4M','E4F'});
Y = [ ...
    % E3M
    mean([E3M_subset.Comp_Hz]) mean([E3M_subset.AMPA_Hz]); ...
    % E3F
    mean([E3F_subset.Comp_Hz]) mean([E3F_subset.AMPA_Hz]); ...
    % E4M
    mean([E4M_subset.Comp_Hz]) mean([E4M_subset.AMPA_Hz]); ...
    % E4F
    mean([E4F_subset.Comp_Hz]) mean([E4F_subset.AMPA_Hz]) ...
    ];
bar(X,Y)
ylabel('Frequency (Hz)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('Compound','AMPA')


%% Plot cumsum between conditions
figure
hold on
plot(mean([E3F_subset.cumsum_event_counts],2),'linewidth',2)
plot(mean([E3M_subset.cumsum_event_counts],2),'linewidth',2)
plot(mean([E4F_subset.cumsum_event_counts],2),'linewidth',2)
plot(mean([E4M_subset.cumsum_event_counts],2),'linewidth',2)
hold off

ylabel('Cumulative event count'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('E3F','E3M','E4F','E4M')


%% Plot baseline changes between conditions

% raw
figure
hold on
plot((mean([E3F_subset.baseline],2))*10e11,'linewidth',2)
plot((mean([E3M_subset.baseline],2))*10e11,'linewidth',2)
plot((mean([E4F_subset.baseline],2))*10e11,'linewidth',2)
plot((mean([E4M_subset.baseline],2))*10e11,'linewidth',2)
hold off
xlabel('Sweep'); ylabel('Amplitude (pA)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('E3F','E3M','E4F','E4M')

% binned

t = 1:180;
t = t';
% E3M
y = [E3M_subset.baseline];
dat = [t y];
[N,edges,bins] = histcounts(dat(:,1),30);
for n = 1:30
bin_means(:,n) = mean(dat(bins==n,:))';
end
E3M_binned_base = mean(bin_means(2:end,:),1);
clear('bin_means')
% E3F
y = [E3F_subset.baseline];
dat = [t y];
[N,edges,bins] = histcounts(dat(:,1),30);
for n = 1:30
bin_means(:,n) = mean(dat(bins==n,:))';
end
E3F_binned_base = mean(bin_means(2:end,:),1);
clear('bin_means')
% E4M
y = [E4M_subset.baseline];
dat = [t y];
[N,edges,bins] = histcounts(dat(:,1),30);
for n = 1:30
bin_means(:,n) = mean(dat(bins==n,:))';
end
E4M_binned_base = mean(bin_means(2:end,:),1);
clear('bin_means')
% E4F
y = [E4F_subset.baseline];
dat = [t y];
[N,edges,bins] = histcounts(dat(:,1),30);
for n = 1:30
bin_means(:,n) = mean(dat(bins==n,:))';
end
E4F_binned_base = mean(bin_means(2:end,:),1);
clear('bin_means')
figure
hold on
plot(E3M_binned_base*10e11,'linewidth',2)
plot(E3F_binned_base*10e11,'linewidth',2)
plot(E4M_binned_base*10e11,'linewidth',2)
plot(E4F_binned_base*10e11,'linewidth',2)
hold off
legend('E3M','E3F','E4M','E4F')
box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Time (mins)')
ylabel('Amplitude (pA)')

% bar
figure;
X = categorical({'E3M','E3F','E4M','E4F'});
X = reordercats(X,{'E3M','E3F','E4M','E4F'});
Y = [ ...
    % E3M
    mean([E3M_subset.early_base]) mean([E3M_subset.late_base]); ...
    % E3F
    mean([E3F_subset.early_base]) mean([E3F_subset.late_base]); ...
    % E4M
    mean([E4M_subset.early_base]) mean([E4M_subset.late_base]); ...
    % E4F
    mean([E4F_subset.early_base]) mean([E4F_subset.late_base]) ...
    ];
err = [ ...
    % E3M
    std([E3M_subset.early_base])/sqrt(length([E3M_subset.early_base])) ...
    std([E3M_subset.late_base])/sqrt(length([E3M_subset.late_base])); ...
    % E3F
    std([E3F_subset.early_base])/sqrt(length([E3F_subset.early_base]))...
    std([E3F_subset.late_base])/sqrt(length([E3F_subset.late_base])); ...
    % E4M
    std([E4M_subset.early_base])/sqrt(length([E4M_subset.early_base]))...
    std([E4M_subset.late_base])/sqrt(length([E4M_subset.late_base])); ...
    % E4F
    std([E4F_subset.early_base])/sqrt(length([E4F_subset.early_base]))...
    std([E4F_subset.late_base])/sqrt(length([E4F_subset.late_base])) ...
    ];
%b = bar(X,Y*10e11,'grouped');

% Grouped error plots
model_series = Y*10e11; 
model_error = err*10e11; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
%Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
ylabel('Amplitude (pA)'); 
box off; set(gca,'linewidth',2); set(gcf,'color','white'), legend('Early Baseline','Late Baseline');
xticklabels({'E3m','E3F','E4M','E4F'})