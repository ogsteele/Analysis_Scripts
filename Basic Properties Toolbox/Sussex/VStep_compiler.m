% VStep_Compiler

% save toggles
save_fig = true;
save_struct = true;
%% Load and compile

% move to correct start point
pathstr = 'D:\Users\Oli Steele\OneDrive\OneDrive - University of Sussex\Ketamine'; 
cd(pathstr)

% list all files of interest
list = dir('**/out.mat'); % ending in .mat
%pat = 'CC_IV_1_full_';
%for i = 1:size(list,1)
%    list(i).isfile = startsWith(list(i).name,pat); % is it one of ours?
%end
%fileFlags = [list.isfile];
%list = list(fileFlags); % save only the ones we're interested in

% create master structure
load(fullfile(list(1).folder,list(1).name)); % load the first file
compiled = out; % rename the first file to facilitate the for loop
% can load file
for i = 2:size(list,1)
    load(fullfile(list(i).folder,list(i).name));
    compiled = [compiled out]; % stitch next file onto the existing file
end

%% Split into condition
E4_control = compiled(strcmp({compiled.Genotype}, 'APOE4') ...
    & strcmp({compiled.Ketamine}, 'No') & strcmp({compiled.Memantine}, 'No'));
E4_ketamine = compiled(strcmp({compiled.Genotype}, 'APOE4') ...
   & strcmp({compiled.Ketamine}, 'Yes') & strcmp({compiled.Memantine}, 'No'));
E3_control = compiled(strcmp({compiled.Genotype}, 'APOE3') ...
    & strcmp({compiled.Ketamine}, 'No') & strcmp({compiled.Memantine}, 'No'));
E3_ketamine = compiled(strcmp({compiled.Genotype}, 'APOE3') ...
    & strcmp({compiled.Ketamine}, 'Yes') & strcmp({compiled.Memantine}, 'No'));
%E3_memantine = compiled(strcmp({compiled.genotype}, 'APOE3') ...
%    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'Yes'));
%E4_memantine = compiled(strcmp({compiled.genotype}, 'APOE4') ...
%    & strcmp({compiled.ketamine}, 'No') & strcmp({compiled.memantine}, 'Yes'));

%% Calculate Maximum Current Densities

ind = [];
% APOE3 Control
for i = 1:size(E3_control,2)
E3C_NaAD(i) = (min(E3_control(i).NaAPeak)*1e12)/E3_control(i).Cap;
E3C_KD(i) = (max(E3_control(i).KMean)*1e12)/E3_control(i).Cap;
E3C_NaID(i) = (min(E3_control(i).NaAPeak)*1e12)/E3_control(i).Cap;
[~,ind(i)] = min(E3_control(i).NaAPeak);
E3C_peak_array(:,i) = E3_control(i).Waves(:,ind(i))/E3_control(i).Cap;
end

% APOE3 Ketamine
for i = 1:size(E3_ketamine,2)
E3K_NaAD(i) = (min(E3_ketamine(i).NaAPeak)*1e12)/E3_ketamine(i).Cap;
E3K_KD(i) = (max(E3_ketamine(i).KMean)*1e12)/E3_ketamine(i).Cap;
E3K_NaID(i) = (min(E3_ketamine(i).NaAPeak)*1e12)/E3_ketamine(i).Cap;
end

ind = [];
% APOE4 Control
for i = 1:size(E4_control,2)
E4C_NaAD(i) = (min(E4_control(i).NaAPeak)*1e12)/E4_control(i).Cap;
E4C_KD(i) = (max(E4_control(i).KMean)*1e12)/E4_control(i).Cap;
E4C_NaID(i) = (min(E4_control(i).NaAPeak)*1e12)/E4_control(i).Cap;
[~,ind(i)] = min(E4_control(i).NaAPeak);
E4C_peak_array(:,i) = E4_control(i).Waves(:,ind(i))/E4_control(i).Cap;
end

% APOE4 Ketamine
for i = 1:size(E4_ketamine,2)
E4K_NaAD(i) = (min(E4_ketamine(i).NaAPeak)*1e12)/E4_ketamine(i).Cap;
E4K_KD(i) = (max(E4_ketamine(i).KMean)*1e12)/E4_ketamine(i).Cap;
E4K_NaID(i) = (min(E4_ketamine(i).NaAPeak)*1e12)/E4_ketamine(i).Cap;
end

%% plot maximum sodium activation current densities
colorscheme = ...
    [0 0 0; ...
    0.6350 0.0780 0.1840];
NaA_ket_fig = figure;
sgtitle('Maximal NaA Current Density:')

    % APOE3 vs APOE4
    labels = {'E3 Control', 'E4 Control'};
    subplot(1,4,1)
    Y = padcat(E3C_NaAD', E4C_NaAD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylim([-400 0])

    % APOE3 with and without ketamine
    labels = {'E3 Control', 'E3 Ketamine'};
    subplot(1,4,2)
    Y = padcat(E3C_NaAD', E3K_NaAD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE3 with and without ketamine')
    ylim([-400 0])

    % APOE4 with and without ketamine
    labels = {'E4 Control', 'E4 Ketamine'};
    subplot(1,4,3)
    Y = padcat(E4C_NaAD', E4K_NaAD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE4 with and without ketamine')
    ylim([-400 0])

    % APOE3/4 with ketamine
    labels = {'E3 Controol', 'E4 Ketamine'};
    subplot(1,4,4)
    Y = padcat(E3C_NaAD', E4K_NaAD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)');
    title('APOE3 vs APOE4 with ketamine')
    ylim([-400 0])

%% plot maximum sodium inactivation current densities
NaI_ket_fig = figure;
sgtitle('Maximal NaI Current Density:')

    % APOE3 vs APOE4
    labels = {'E3 Control', 'E4 Control'};
    subplot(1,4,1)
    Y = padcat(E3C_NaID', E4C_NaID');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)');
    title('APOE3 vs APOE4')
    ylim([-400 0])

    % APOE3 with and without ketamine
    labels = {'E3 Control', 'E3 Ketamine'};
    subplot(1,4,2)
    Y = padcat(E3C_NaID', E3K_NaID');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE3 with and without ketamine')
    ylim([-400 0])

    % APOE4 with and without ketamine
    labels = {'E4 Control', 'E4 Ketamine'};
    subplot(1,4,3)
    Y = padcat(E4C_NaID', E4K_NaID');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE4 with and without ketamine')
    ylim([-400 0])

    % APOE3/4 with ketamine
    labels = {'E3 Control', 'E4 Ketamine'};
    subplot(1,4,4)
    Y = padcat(E3C_NaID', E4K_NaID');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)');
    title('APOE3 vs APOE4 with ketamine')
    ylim([-400 0])

%% plot maximum potassium current densities
K_ket_fig = figure;
sgtitle('Maximal K Current Density:')

    % APOE3 vs APOE4
    labels = {'E3 Control', 'E4 Control'};
    subplot(1,4,1)
    Y = padcat(E3C_KD', E4C_KD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)');
    title('APOE3 vs APOE4')
    ylim([0 200])

    % APOE3 with and without ketamine
    labels = {'E3 Control', 'E3 Ketamine'};
    subplot(1,4,2)
    Y = padcat(E3C_KD', E3K_KD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE3 with and without ketamine')
    ylim([0 200])

    % APOE4 with and without ketamine
    labels = {'E4 Control', 'E4 Ketamine'};
    subplot(1,4,3)
    Y = padcat(E4C_KD', E4K_KD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)')
    title('APOE4 with and without ketamine')
    ylim([0 200])

    % APOE3/4 with ketamine
    labels = {'E3 Control', 'E4 Ketamine'};
    subplot(1,4,4)
    Y = padcat(E3C_KD', E4K_KD');
    b = boxplot(Y, 'labels', labels);
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    set(b(7,:),'Visible','off') % make the outlier points invisible
    x=repmat(1:4,length(Y),1);
    hold on
    % plot individual values
    for i = 1:size(Y,2)
        scatter(x(:,i),Y(:,i),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15,'MarkerFacecolor',colorscheme(i,:));
    end
    ylabel('Current Density (pA/pF)');
    title('APOE3 vs APOE4 with ketamine')
    ylim([0 200])

%% Current voltage relationships


%% save figures and output

% navigate to output directory
cd('VStep Analysis')

% structure
vals.E3C_NaAD = E3C_NaAD;
vals.E4C_NaAD = E4C_NaAD;
vals.E3C_KD = E3C_KD;
vals.E4C_KD = E4C_KD;
vals.E3C_peak_array = E3C_peak_array;
vals.E4C_peak_array = E4C_peak_array;

if save_struct == true
save('compiled.mat', 'compiled', '-v7.3')
save('vals.mat', 'vals', '-v7.3')
else
end

% figures
if save_fig == true
%savefig(NaA_fig, 'NaA_fig')
savefig(NaA_ket_fig, 'NaA_ket_fig')
%savefig(K_fig,'K_fig')
savefig(K_ket_fig,'K_ket_fig')
%savefig(NaI_fig,'NaI_fig')
savefig(NaI_ket_fig, 'NaI_ket_fig')
else 
end
