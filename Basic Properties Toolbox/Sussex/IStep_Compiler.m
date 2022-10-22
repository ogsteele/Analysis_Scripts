% IStep_Comp

%% Fig Toggle
mem_tog = 1;
exc_tog = 1;
iei_tog = 0;
rh_tog = 1;

%% Code

% move to correct start point
pathstr = 'E:\Experiments\Ketamine mEPSC\04_processing'; 
cd(pathstr)

% list all files of interest
list = dir('**\*.mat'); % ending in .mat
pat = 'CC_IV_1_full_';
for i = 1:size(list,1)
    list(i).isfile = startsWith(list(i).name,pat); % is it one of ours?
end
fileFlags = [list.isfile];
list = list(fileFlags); % save only the ones we're interested in

% create master structure
load(fullfile(list(1).folder,list(1).name)); % load the first file
compiled = output; % rename the first file to facilitate the for loop
% can load file
for i = 2:size(list,1)
    load(fullfile(list(i).folder,list(i).name));
    compiled = [compiled output]; % stitch next file onto the existing file
end

% correct for input mismatch
x = zeros(31,1); % empty array of cells
for i = 1:size(compiled,2)
    if size(compiled(i).numSpikes,1) == 30
        spikesHolder = x; % assign holder variable the empty array
        spikesHolder(2:end) = compiled(i).numSpikes;
        compiled(i).numSpikes = spikesHolder;
    else
    end
end

%% Penn IEI Setup

if iei_tog == 1
Waves = compiled(21).waves;
% preallocate pks, locs, w, p, numSpikes and mute error here
warning('off','signal:findpeaks:largeMinPeakHeight');

locs = cell(size(Waves,2),1);
difflocs = cell(size(Waves,2),1); % difference between loc in dp
normlocs = cell(size(Waves,2),1); % difference between loc in dp
nlocs = cell(size(Waves,2),1); % interval index
% 
% Y = padcat([E3_control.(features{j})]',[E3_ketamine.(features{j})]',[E4_control.(features{j})]',[E4_ketamine.(features{j})]');
% 
for s = 1:size(compiled,2)
% findpeaks to determine number and location of AP's per wave
    % create blank vals for nested for loop
    compiled(s).interval_index = [];
    compiled(s).normalised_interval = [];
    for i = 1:size(compiled(s).waves,2)
        % pks = value of peak
        % locs = index of peak
        % w = width
        % p = prominence
        [~,compiled(s).locs{i},~,p] = findpeaks(compiled(s).waves(:,i),...
            'MinPeakHeight',0,...
            'MaxPeakWidth',0.02*10^4,...
            'MinPeakProminence',0.01,...
            'MinPeakDistance',0.02*10^4,...
            'WidthReference','halfheight');
        compiled(s).difflocs{i} = diff(compiled(s).locs{i}); % difference between loc in dp
        compiled(s).normlocs{i} = normalize(compiled(s).difflocs{i},'scale','first'); % difflocs norm to first value
        compiled(s).nlocs{i} = 1:size(compiled(s).difflocs{i}); % interval index
        compiled(s).interval_index = [compiled(s).interval_index,compiled(s).nlocs{i}];
        compiled(s).normalised_interval = [compiled(s).normalised_interval,compiled(s).normlocs{i}'];
    end
end
else
    disp('iei_tog == 0')
end

%% Split into condition
E4_control = compiled(strcmp({compiled.genotype}, 'APOE4') ...
    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'No'));
E4_ketamine = compiled(strcmp({compiled.genotype}, 'APOE4') ...
    & strcmp({compiled.ketamine}, 'Yes') & strcmp({compiled.memantine}, 'No'));
E3_control = compiled(strcmp({compiled.genotype}, 'APOE3') ...
    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'No'));
E3_ketamine = compiled(strcmp({compiled.genotype}, 'APOE3') ...
    & strcmp({compiled.ketamine}, 'Yes') & strcmp({compiled.memantine}, 'No'));
E3_memantine = compiled(strcmp({compiled.genotype}, 'APOE3') ...
    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'Yes'));
E4_memantine = compiled(strcmp({compiled.genotype}, 'APOE4') ...
    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'Yes'));

%% plotting variables
colorscheme = ...
    [0 0 0; ...
    0 0 0; ...
    0.6350 0.0780 0.1840; ...
    0.6350 0.0780 0.1840];
labels = ...
    {'E3 Control', 'E3 Ketamine', 'E4 Control', 'E4 Ketamine'};
features = ...
    {'Rh','peak','afterhyp','amp','thresh','half','rise','fall','IR'};
ylabels = ...
    {'pA','mV','mV','mV','mV','ms','mV/ms','mV/ms','MOhm'};

%% summary data plots
if mem_tog == 1
figure;
for j = 1:size(features,2)
    subplot(3,3,j)
    Y = padcat([E3_control.(features{j})]',[E3_ketamine.(features{j})]',[E4_control.(features{j})]',[E4_ketamine.(features{j})]');
    b = boxplot(Y,...
            'labels',labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    hold off; ylabel(ylabels{j}); title(features{j})
end
else
    disp('mem_tog == 0')
end

%% membrane excitability plots
if exc_tog == 1
aveSpikes = mean([E3_control.numSpikes],2); % mean number of spikes per I
conditionList = ...
    {'E3_control', 'E3_ketamine', 'E4_control', 'E4_ketamine'};

for i = 1:size(E3_control,2)
    E3C_peakSpikes(i) = max(E3_control(i).numSpikes); % maximum number of APs fired
    E3C_totalSpikes(i) = sum(E3_control(i).numSpikes); % total number of APs fired
end
E3C_aveSpikes = mean([E3_control.numSpikes],2); % mean number of spikes per I
x = [E3_control.numSpikes];

for i = 1:size(x,1)
    E3C_aveSpikesSTD(i) = std(x(i,:));
    E3C_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for  i = size(E3_ketamine,2)
    E3K_peakSpikes(i) = max(E3_ketamine(i).numSpikes); % maximum number of APs fired
    E3K_totalSpikes(i) = sum(E3_ketamine(i).numSpikes); % total number of APs fired
end
E3K_aveSpikes = mean([E3_ketamine.numSpikes],2); % mean number of spikes per I
x = [E3_ketamine.numSpikes];

for i = 1:size(x,1)
    E3K_aveSpikesSTD(i) = std(x(i,:));
    E3K_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for i = 1:size(E4_control,2)
    E4C_peakSpikes(i) = max(E4_control(i).numSpikes); % maximum number of APs fired
    E3C_totalSpikes(i) = sum(E4_control(i).numSpikes); % total number of APs fired
end
E4C_aveSpikes = mean([E4_control.numSpikes],2); % mean number of spikes per I
x = [E4_control.numSpikes];
for i = 1:size(x,1)
    E4C_aveSpikesSTD(i) = std(x(i,:));
E4C_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for i = 1:size(E4_ketamine,2)
    E4K_peakSpikes(i) = max(E4_ketamine(i).numSpikes); % maximum number of APs fired
    E4K_totalSpikes(i) = sum(E4_ketamine(i).numSpikes); % total number of APs fired
end
E4K_aveSpikes = mean([E4_ketamine.numSpikes],2); % mean number of spikes per I
x = [E4_ketamine.numSpikes];
for i = 1:size(x,1)
    E4K_aveSpikesSTD(i) = std(x(i,:));
E4K_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

pA = [-200;-180;-160;-140;-120;-100;-80;-60;-40;-20;0;20;40;60;80;100;120;140;160;180;200;220;240;260;280;300;320;340;360;380;400];

figure;
subplot(1,4,1); plot(pA,E3C_aveSpikes,'color',[0 0.4470 0.7410]); 
hold on; plot(pA,E4C_aveSpikes,'color',[0.8500 0.3250 0.0980]);
errorbar(pA,E3C_aveSpikes,E3C_aveSpikesSEM,'HandleVisibility','off','color',[0 0.4470 0.7410])
errorbar(pA,E4C_aveSpikes,E4C_aveSpikesSEM,'HandleVisibility','off','color',[0.8500 0.3250 0.0980])
legend('E3','E4','linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('APOE3 vs APOE4')
ylim([0 14])

subplot(1,4,2); plot(pA,E3C_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E3K_aveSpikes,'color',[0.6350 0.0780 0.1840]);
errorbar(pA,E3C_aveSpikes,E3C_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E3K_aveSpikes,E3K_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend(labels(1:2),'linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E3 +/- Ketamine')
ylim([0 14])

subplot(1,4,3); plot(pA,E4C_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E4K_aveSpikes,'color',[0.6350 0.0780 0.1840]); 
errorbar(pA,E4C_aveSpikes,E4C_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E4K_aveSpikes,E4K_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend(labels(3:4),'linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E4 +/- Ketamine')
ylim([0 14])

subplot(1,4,4); plot(pA,E3K_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E4K_aveSpikes,'color',[0.6350 0.0780 0.1840]); 
errorbar(pA,E3K_aveSpikes,E3K_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E4K_aveSpikes,E4K_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend('E3 Ketamine','E4 Ketamine','linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E3 vs E4 + Ketamine')
ylim([0 14])
else 
    disp('exc_tog == 0');
end

%% membrane excitability plots (memantine)
if exc_tog == 1
memlabels =  {'E3 Control', 'E3 Memantine', 'E4 Control', 'E4 Memantine'};
aveSpikes = mean([E3_control.numSpikes],2); % mean number of spikes per I
conditionList = ...
    {'E3_control', 'E3_memantine', 'E4_control', 'E4_memantine'};

for i = 1:size(E3_control,2)
    E3C_peakSpikes(i) = max(E3_control(i).numSpikes); % maximum number of APs fired
    E3C_totalSpikes(i) = sum(E3_control(i).numSpikes); % total number of APs fired
end
E3C_aveSpikes = mean([E3_control.numSpikes],2); % mean number of spikes per I
x = [E3_control.numSpikes];
for i = 1:size(x,1)
    E3C_aveSpikesSTD(i) = std(x(i,:));
    E3C_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for  i = size(E3_memantine,2)
    E3M_peakSpikes(i) = max(E3_memantine(i).numSpikes); % maximum number of APs fired
    E3M_totalSpikes(i) = sum(E3_memantine(i).numSpikes); % total number of APs fired
end
E3M_aveSpikes = mean([E3_memantine.numSpikes],2); % mean number of spikes per I
x = [E3_memantine.numSpikes];
for i = 1:size(x,1)
    E3M_aveSpikesSTD(i) = std(x(i,:));
    E3M_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for i = 1:size(E4_control,2)
    E4C_peakSpikes(i) = max(E4_control(i).numSpikes); % maximum number of APs fired
    E3C_totalSpikes(i) = sum(E4_control(i).numSpikes); % total number of APs fired
end
E4C_aveSpikes = mean([E4_control.numSpikes],2); % mean number of spikes per I
x = [E4_control.numSpikes];
for i = 1:size(x,1)
    E4C_aveSpikesSTD(i) = std(x(i,:));
E4C_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

for i = 1:size(E4_memantine,2)
    E4M_peakSpikes(i) = max(E4_memantine(i).numSpikes); % maximum number of APs fired
    E4M_totalSpikes(i) = sum(E4_memantine(i).numSpikes); % total number of APs fired
end
E4M_aveSpikes = mean([E4_memantine.numSpikes],2); % mean number of spikes per I
x = [E4_memantine.numSpikes];
for i = 1:size(x,1)
    E4M_aveSpikesSTD(i) = std(x(i,:));
    E4M_aveSpikesSEM(i) = std(x(i,:))/sqrt(length(x(i,:)));
end

pA = [-200;-180;-160;-140;-120;-100;-80;-60;-40;-20;0;20;40;60;80;100;120;140;160;180;200;220;240;260;280;300;320;340;360;380;400];

figure;
subplot(1,4,1); plot(pA,E3C_aveSpikes,'color',[0 0.4470 0.7410]); 
hold on; plot(pA,E4C_aveSpikes,'color',[0.8500 0.3250 0.0980]);
errorbar(pA,E3C_aveSpikes,E3C_aveSpikesSEM,'HandleVisibility','off','color',[0 0.4470 0.7410])
errorbar(pA,E4C_aveSpikes,E4C_aveSpikesSEM,'HandleVisibility','off','color',[0.8500 0.3250 0.0980])
legend('E3','E4','linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('APOE3 vs APOE4')
ylim([0 14])

subplot(1,4,2); plot(pA,E3C_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E3M_aveSpikes,'color',[0.6350 0.0780 0.1840]);
errorbar(pA,E3C_aveSpikes,E3C_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E3M_aveSpikes,E3M_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend(memlabels(1:2),'linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E3 +/- Memantine')
ylim([0 14])

subplot(1,4,3); plot(pA,E4C_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E4M_aveSpikes,'color',[0.6350 0.0780 0.1840]); 
errorbar(pA,E4C_aveSpikes,E4C_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E4M_aveSpikes,E4M_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend(memlabels(3:4),'linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E4 +/- Memantine')
ylim([0 14])

subplot(1,4,4); plot(pA,E3M_aveSpikes,'color',[0 0 0]); 
hold on; plot(pA,E4M_aveSpikes,'color',[0.6350 0.0780 0.1840]); 
errorbar(pA,E3M_aveSpikes,E3M_aveSpikesSEM,'HandleVisibility','off','color',[0 0 0])
errorbar(pA,E4M_aveSpikes,E4M_aveSpikesSEM,'HandleVisibility','off','color',[0.6350 0.0780 0.1840])
legend('E3 Memantine','E4 Memantine','linewidth',1,'location','northwest'); box off; set(gca,'linewidth',2); set(gcf,'color','white')
xlabel('Current Injected (pA)'); ylabel('Number of Action Potentials'); title('E3 vs E4 + Memantine')
ylim([0 14])
else 
    disp('exc_tog == 0');
end

%% Penn IEI Plots
if iei_tog == 1
% E3 Events
% create blank vals for nested for loop
interval_index = [];
normalised_interval = [];
for s = 1:size(E3_control,2) 
    interval_index = [interval_index,E3_control(s).interval_index];
    normalised_interval = [normalised_interval,E3_control(s).normalised_interval];
end
E3_interval_index = interval_index;
E3_normalised_interval = normalised_interval;

% E3K Events
% create blank vals for nested for loop
interval_index = [];
normalised_interval = [];
for s = 1:size(E3_ketamine,2) 
    interval_index = [interval_index,E3_ketamine(s).interval_index];
    normalised_interval = [normalised_interval,E3_ketamine(s).normalised_interval];
end
E3K_interval_index = interval_index;
E3K_normalised_interval = normalised_interval;
    
% E4 events
% create blank vals for nested for loop
interval_index = [];
normalised_interval = [];
for s = 1:size(E4_control,2) 
    interval_index = [interval_index,E4_control(s).interval_index];
    normalised_interval = [normalised_interval,E4_control(s).normalised_interval];
end
E4_interval_index = interval_index;
E4_normalised_interval = normalised_interval;


% E4K events
% create blank vals for nested for loop
interval_index = [];
normalised_interval = [];
for s = 1:size(E4_ketamine,2)
    interval_index = [interval_index,E4_ketamine(s).interval_index];
    normalised_interval = [normalised_interval,E4_ketamine(s).normalised_interval];
end
E4K_interval_index = interval_index;
E4K_normalised_interval = normalised_interval;

% coefficients, fits and errors
[E3p,E3S] = polyfit(E3_interval_index,E3_normalised_interval,1); % coefficients of first degree polynomial
[E3y_fit,E3delta] = polyval(E3p,E3_interval_index,E3S); % fit the polynomial and error

[E3Kp,E3KS] = polyfit(E3K_interval_index,E3K_normalised_interval,1); % coefficients of first degree polynomial
[E3Ky_fit,E3Kdelta] = polyval(E3Kp,E3K_interval_index,E3KS); % fit the polynomial and error

[E4p,E4S] = polyfit(E4_interval_index,E4_normalised_interval,1); % coefficients of first degree polynomial
[E4y_fit,E4delta] = polyval(E4p,E4_interval_index,E4S); % fit the polynomial and error

[E4Kp,E4KS] = polyfit(E4K_interval_index,E4K_normalised_interval,1); % coefficients of first degree polynomial
[E4Ky_fit,E4Kdelta] = polyval(E4Kp,E4K_interval_index,E4KS); % fit the polynomial and error


figure;subplot(1,4,1)
scatter(E3_interval_index,E3_normalised_interval,'bo','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
hold on
plot(E3_interval_index,E3y_fit,'b-','linewidth',3)
%plot(E3_interval_index,E3y_fit+2*E3delta,'b--',E3_interval_index,E3y_fit-2*E3delta,'b--')
scatter(E4_interval_index,E4_normalised_interval,'ro','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
plot(E4_interval_index,E4y_fit,'r-','linewidth',3)
%plot(E4_interval_index,E4y_fit+2*E4delta,'r--',E4_interval_index,E4y_fit-2*E4delta,'r--')
legend('E3 Data','Linear Fit','E4 Data','Linear Fit')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('AP Interval #'); ylabel('Normalised Inter Event Interval')
title('APOE3 vs APOE4')
ylim([0 15]); xlim([0 15]);

subplot(1,4,2)
scatter(E3_interval_index,E3_normalised_interval,'bo','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
hold on
plot(E3_interval_index,E3y_fit,'b-','linewidth',3)
%plot(E3_interval_index,E3y_fit+2*E3delta,'b--',E3_interval_index,E3y_fit-2*E3delta,'b--')
scatter(E3K_interval_index,E3K_normalised_interval,'ro','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
plot(E3K_interval_index,E3Ky_fit,'r-','linewidth',3)
%plot(E3K_interval_index,E3Ky_fit+2*E3Kdelta,'r--',E3K_interval_index,E3Ky_fit-2*E3Kdelta,'r--')
legend('E3 Data','Linear Fit','E3K Data','Linear Fit')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('AP Interval #'); ylabel('Normalised Inter Event Interval')
title('APOE3 +/- Ketamine')
ylim([0 15]); xlim([0 15]);

subplot(1,4,3)
scatter(E4_interval_index,E4_normalised_interval,'bo','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
hold on
plot(E4_interval_index,E4y_fit,'b-','linewidth',3)
%plot(E4_interval_index,E4y_fit+2*E4delta,'b--',E4_interval_index,E4y_fit-2*E4delta,'b--')
scatter(E4K_interval_index,E4K_normalised_interval,'ro','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
plot(E4K_interval_index,E4Ky_fit,'r-','linewidth',3)
%plot(E4K_interval_index,E4Ky_fit+2*E4Kdelta,'r--',E4K_interval_index,E4Ky_fit-2*E4Kdelta,'r--')
legend('E4 Data','Linear Fit','E4K Data','Linear Fit')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('AP Interval #'); ylabel('Normalised Inter Event Interval')
title('APOE4 +/- Ketamine')
ylim([0 15]); xlim([0 15]);

subplot(1,4,4)
scatter(E3_interval_index,E3_normalised_interval,'bo','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
hold on
plot(E3_interval_index,E3y_fit,'b-','linewidth',3)
%plot(E3_interval_index,E3y_fit+2*E3delta,'b--',E3_interval_index,E3y_fit-2*E3delta,'b--')
scatter(E4K_interval_index,E4K_normalised_interval,'ro','filled','MarkerFaceAlpha',0.1,'jitter','on', 'jitterAmount',0.5)
plot(E4K_interval_index,E4Ky_fit,'r-','linewidth',3)
%plot(E4K_interval_index,E4Ky_fit+2*E4Kdelta,'r--',E4K_interval_index,E4Ky_fit-2*E4Kdelta,'r--')
legend('E3 Data','Linear Fit','E4K Data','Linear Fit')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('AP Interval #'); ylabel('Normalised Inter Event Interval')
title('E3 vs E4 + Ketamine')
ylim([0 15]); xlim([0 15]);

else
    disp('iei_tog == 0')
end


%% Alt_PennIEI Plot
% plot the AP # against the inter event interval

% first get the very first IEI from all conditions
wave_ind = zeros(size(compiled,2),6);
pA_step = zeros(size(compiled,2),6);
init_interval = zeros(size(compiled,2),6);
pA = [-200;-180;-160;-140;-120;-100;-80;-60;-40;-20;0;20;40;60;80;100;120;140;160;180;200;220;240;260;280;300;320;340;360;380;400];
for j = 1:6
    for i = 1:size(compiled,2)
        temp = [];
        wave_ind(i,j) = (find(~cellfun(@isempty,compiled(i).difflocs),1))+(j-1); % first wave with 2x AP
        pA_step(i,j) = pA(wave_ind(i,j)); % current step for the first wave with two APs
        temp = compiled(i).difflocs{wave_ind(i,j)}; % holding var of the first val
        if isempty(temp) == 1
            init_interval(i,j) = NaN;
        else
            init_interval(i,j) = temp(1); % the initial interval across all conditions
        end
    end
end

%% Rheobase Plot
if rh_tog == 1
Y = padcat([E3_control.Rh]',[E3_ketamine.Rh]',[E4_control.Rh]',[E4_ketamine.Rh]');
b = boxplot(Y,...
        'labels',labels);
box off; set(gca,'linewidth',2); set(gcf,'color','white');
set(b(7,:),'Visible','off') % make the outlier points invisible
x=repmat(1:4,length(Y),1);
hold on
% plot individual values
for i = 1:size(Y,2)
    scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
end
xlabel('Condition'); ylabel('Rheobase (pA)');
else
    disp('rh_tog == 0')
end

%% Phase profiles

[E3C_mv,E3C_dvdt] = posthocPhaseprofile(E3_control);
[E4C_mv,E4C_dvdt] = posthocPhaseprofile(E4_control);
[E4K_mv,E4K_dvdt] = posthocPhaseprofile(E4_ketamine);
[E4M_mv,E4M_dvdt] = posthocPhaseprofile(E4_memantine);

figure; plot(mean(E3C_mv,2),mean(E3C_dvdt,2),...
    'color','black','linewidth',2,...
    'DisplayName','E3C')
hold on; plot(mean(E4C_mv,2),mean(E4C_dvdt,2),...
    'color','red','linewidth',2,...
    'DisplayName','E4C')
box off; set(gca,'linewidth',2); set(gcf,'color','white')
title('Control'); xlabel('V (mV'); ylabel('dvdt (mV/ms'); legend('LineWidth',1)

figure; plot(mean(E3C_mv,2),mean(E3C_dvdt,2),...
    'color','black','linewidth',2,...
    'DisplayName','E3C')
hold on; plot(mean(E4K_mv,2),mean(E4K_dvdt,2),...
    'color','red','linewidth',2,...
    'DisplayName','E4K')
box off; set(gca,'linewidth',2); set(gcf,'color','white')
title('Ketamine'); xlabel('V (mV'); ylabel('dvdt (mV/ms'); legend('LineWidth',1)

figure; plot(mean(E3C_mv,2),mean(E3C_dvdt,2),...
    'color','black','linewidth',2,...
    'DisplayName','E3C')
hold on; plot(mean(E4M_mv,2),mean(E4M_dvdt,2),...
    'color','red','linewidth',2,...
    'DisplayName','E4M')
box off; set(gca,'linewidth',2); set(gcf,'color','white')
title('Memantine'); xlabel('V (mV'); ylabel('dvdt (mV/ms'); legend('LineWidth',1)

%% Save Variables
save('compiled.mat', 'compiled', '-v7.3')
save conditions.mat E3_ketamine E3_control E4_ketamine E4_control E3_memantine E4_memantine
clearvars -except compiled E3_ketamine E3_control E4_ketamine E4_control E3_memantine E4_memantine
load compiled.mat

% ------
function [mV,dvdt] = posthocPhaseprofile(condition)

% https://www.nature.com/articles/nrn2148
% https://www.cell.com/biophysj/pdf/S0006-3495(63)86827-7.pdf
% Tempted not to be too distracted by the phase plane plots 
% however as Bean et al suggest they're less informative when
% recorded from cell bodies due the confounding presence of the
% dendritic tree / axon. Not to mention these are also action 
% potentials after current injection, which further 
% complicates it

% post hoc phase profile of the AP script
% to be included in the IStep_compiler graph

% set the detection zones
S.xdiff = condition(1).ephysIO.xdiff;
detection = [2*1e4,6*1e4];
detStart = round(detection(1)/S.xdiff);

% set the winSize_ms and LPF_Hz variables
winSize_ms = 25;
winSize_ms = (winSize_ms*1e-3)/S.xdiff; % (ms to s)to data points
LPF_Hz = 330;

% figure; % figure when debugging
% construct for loop to run process
for i = 1:size(condition,2)
    % determine the wavenumber of first AP
    locs = condition(i).locs;
        C = locs;
            idx = ~cellfun('isempty',C);
            outs = zeros(size(C));
            outs(idx) = cellfun(@(v)v(1),C(idx));
            for lp = 1:size(outs,1)
                if outs(lp) > 0
                outs(lp) = outs(lp) + detStart; % account for the detection zone
                end
            end
            logicalIndexes =  outs < 27500 & outs > 1;
    wavenum_first(i) = (size(locs,2)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
   
    % detemine AP window for wave of interest
    warning('off','MATLAB:colon:nonIntegerIndex')
    mV(:,i) = condition(i).waves(outs(wavenum_first(i))-(winSize_ms/2):outs(wavenum_first(i))+(winSize_ms/2),wavenum_first(i))*1000;
    % calcuate dv/dt
    dvdt(:,i) = gradient(mV(:,i));
    warning('on','MATLAB:colon:nonIntegerIndex')
    %hold on; plot(mV(:,i)); legend() % plot when debugging
end
end
