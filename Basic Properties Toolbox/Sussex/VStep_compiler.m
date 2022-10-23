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
E3_memantine = compiled(strcmp({compiled.Genotype}, 'APOE3') ...
    & strcmp({compiled.Ketamine}, 'No') & strcmp({compiled.Memantine}, 'Yes'));
E4_memantine = compiled(strcmp({compiled.Genotype}, 'APOE4') ...
    & strcmp({compiled.Ketamine}, 'No') & strcmp({compiled.Memantine}, 'Yes'));

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
[~,ind(i)] = min(E3_ketamine(i).NaAPeak);
E3K_peak_array(:,i) = E3_ketamine(i).Waves(:,ind(i))/E3_ketamine(i).Cap;
end

% APOE3 Memantine
for i = 1:size(E3_memantine,2)
E3M_NaAD(i) = (min(E3_memantine(i).NaAPeak)*1e12)/E3_memantine(i).Cap;
E3M_KD(i) = (max(E3_memantine(i).KMean)*1e12)/E3_memantine(i).Cap;
E3M_NaID(i) = (min(E3_memantine(i).NaAPeak)*1e12)/E3_memantine(i).Cap;
[~,ind(i)] = min(E3_memantine(i).NaAPeak);
E3M_peak_array(:,i) = E3_memantine(i).Waves(:,ind(i))/E3_memantine(i).Cap;
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

% APOE4 Ketaminee
for i = 1:size(E4_ketamine,2)
E4K_NaAD(i) = (min(E4_ketamine(i).NaAPeak)*1e12)/E4_ketamine(i).Cap;
E4K_KD(i) = (max(E4_ketamine(i).KMean)*1e12)/E4_ketamine(i).Cap;
E4K_NaID(i) = (min(E4_ketamine(i).NaAPeak)*1e12)/E4_ketamine(i).Cap;
[~,ind(i)] = min(E4_ketamine(i).NaAPeak);
E4K_peak_array(:,i) = E4_ketamine(i).Waves(:,ind(i))/E4_ketamine(i).Cap;
end

% APOE4 Memantine
for i = 1:size(E4_memantine,2)
E4M_NaAD(i) = (min(E4_memantine(i).NaAPeak)*1e12)/E4_memantine(i).Cap;
E4M_KD(i) = (max(E4_memantine(i).KMean)*1e12)/E4_memantine(i).Cap;
E4M_NaID(i) = (min(E4_memantine(i).NaAPeak)*1e12)/E4_memantine(i).Cap;
[~,ind(i)] = min(E4_memantine(i).NaAPeak);
E4M_peak_array(:,i) = E4_memantine(i).Waves(:,ind(i))/E4_memantine(i).Cap;
end

%% plot maximum sodium activation current densities +/- ketamine
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

%% plot maximum sodium activation current densities +/- memantine
colorscheme = ...
    [0 0 0; ...
    0.6350 0.0780 0.1840];
NaA_mem_fig = figure;
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

    % APOE3 with and without memantine
    labels = {'E3 Control', 'E3 Memantine'};
    subplot(1,4,2)
    Y = padcat(E3C_NaAD', E3M_NaAD');
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
    title('APOE3 with and without memantine')
    ylim([-400 0])

    % APOE4 with and without Memantine
    labels = {'E4 Control', 'E4 Memantine'};
    subplot(1,4,3)
    Y = padcat(E4C_NaAD', E4M_NaAD');
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
    title('APOE4 with and without memantine')
    ylim([-400 0])

    % APOE3/4 with memantine
    labels = {'E3 Controol', 'E4 Ketamine'};
    subplot(1,4,4)
    Y = padcat(E3C_NaAD', E4M_NaAD');
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
    title('APOE3 vs APOE4 with memantine')
    ylim([-400 0])

%% plot maximum sodium inactivation current densities +/- ketamine
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

%% plot maximum sodium inactivation current densities +/- memantine
NaI_mem_fig = figure;
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

    % APOE3 with and without memantine
    labels = {'E3 Control', 'E3 Memantine'};
    subplot(1,4,2)
    Y = padcat(E3C_NaID', E3M_NaID');
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
    title('APOE3 with and without Memantine')
    ylim([-400 0])

    % APOE4 with and without Memantine
    labels = {'E4 Control', 'E4 Memantine'};
    subplot(1,4,3)
    Y = padcat(E4C_NaID', E4M_NaID');
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
    title('APOE4 with and without Memantine')
    ylim([-400 0])

    % APOE3/4 with Memantine
    labels = {'E3 Control', 'E4 Memantine'};
    subplot(1,4,4)
    Y = padcat(E3C_NaID', E4M_NaID');
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
    title('APOE3 vs APOE4 with Memantine')
    ylim([-400 0])

%% plot maximum potassium current densities +/- ketamine
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

%% plot maximum potassium current densities +/- Memantine
K_mem_fig = figure;
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

    % APOE3 with and without Memantine
    labels = {'E3 Control', 'E3 Memantine'};
    subplot(1,4,2)
    Y = padcat(E3C_KD', E3M_KD');
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
    title('APOE3 with and without Memantine')
    ylim([0 200])

    % APOE4 with and without Memantine
    labels = {'E4 Control', 'E4 Memantine'};
    subplot(1,4,3)
    Y = padcat(E4C_KD', E4M_KD');
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
    title('APOE4 with and without Memantine')
    ylim([0 200])

    % APOE3/4 with ketamine
    labels = {'E3 Control', 'E4 Memantine'};
    subplot(1,4,4)
    Y = padcat(E3C_KD', E4M_KD');
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
    title('APOE3 vs APOE4 with Memantine')
    ylim([0 200])

%% Peak NaA waves for trace overlays

% baseline subtract each of these arrays independently
% E3C
baseline = mean(E3C_peak_array(1:500,:));
E3C_peak_array = E3C_peak_array - baseline;
% E3K
baseline = mean(E3K_peak_array(1:500,:));
E3K_peak_array = E3K_peak_array - baseline;
% E3M
baseline = mean(E3M_peak_array(1:500,:));
E3M_peak_array = E3M_peak_array - baseline;
% E4C
baseline = mean(E4C_peak_array(1:500,:));
E4C_peak_array = E4C_peak_array - baseline;
% E4K
baseline = mean(E4K_peak_array(1:500,:));
E4K_peak_array = E4K_peak_array - baseline;
% E4M
baseline = mean(E4M_peak_array(1:500,:));
E4M_peak_array = E4M_peak_array - baseline;

% peak align traces
detZone = [1000,4000];
alignedE3C = peakAlign(E3C_peak_array,detZone);
alignedE4C = peakAlign(E4C_peak_array,detZone);
alignedE3K = peakAlign(E3K_peak_array,detZone);
alignedE4K = peakAlign(E4K_peak_array,detZone);
alignedE3M = peakAlign(E3M_peak_array,detZone);
alignedE4M = peakAlign(E4M_peak_array,detZone);

% calculate averages of each wave, and pass through a filter
meanE3C = mean(alignedE3C,2);
meanE4C = mean(alignedE4C,2);
meanE3K = mean(alignedE3K,2);
meanE4K = mean(alignedE4K,2);
meanE3M = mean(alignedE3M,2);
meanE4M = mean(alignedE4M,2);

% create exemplary overlay plots of the regions of interest
NaC_overlay = figure; plot(meanE3C,'color','black','DisplayName','E3C')
hold on; plot(meanE4C,'color','red','DisplayName','E4C')
legend

NaK_overlay = figure; plot(meanE3C,'color','black','DisplayName','E3C')
hold on; plot(meanE4K,'color','red','DisplayName','E4K')
legend

NaM_overlay = figure; plot(meanE3C,'color','black','DisplayName','E3C')
hold on; plot(meanE4M,'color','red','DisplayName','E4M')
legend



%% Current voltage relationships
% note, not sure what these are truly supposed to be

IDplots(E3_control,'APOE3',E4_control,'APOE4')
IDplots(E3_control,'APOE3',E4_ketamine,'APOE4 + Ket')
IDplots(E3_control,'APOE3',E4_memantine,'APOE4+ Mem')

%     convert to pA and divide by pF to calculate the current density in
%     in pA/pF
%     NaK_ID_fig = figure;
%     NaID = (NaApeakval*1e12)/capacitance;
%     KID = (Kmeanval*1e12)/capacitance;
%     plot(Vm,NaID,'-ob')
%     hold on; plot(Vm,KID,'-or')
%     ax = gca;
%     ax.XAxisLocation = 'origin'; xlabel('Membrane Potential (mV)')
%     ax.YAxisLocation = 'origin'; ylabel('Current Density (pA/pF)')
%     box off; set(gcf,'color','white'); set(gca,'linewidth',2)
%     legend('Sodium Activation','Potassium','location','southwest','linewidth',1);
%     sgtitle('Na/K Current Density')

%% save figures and output

% navigate to output directory
cd('VStep Analysis')

% structure
vals.E3C_NaAD = E3C_NaAD;
vals.E4C_NaAD = E4C_NaAD;
vals.E3K_NaAD = E3K_NaAD;
vals.E4K_NaAD = E4K_NaAD;
vals.E3M_NaAD = E3M_NaAD;
vals.E4M_NaAD = E4M_NaAD;
vals.E3C_KD = E3C_KD;
vals.E4C_KD = E4C_KD;
vals.E3K_KD = E3K_KD;
vals.E4K_KD = E4K_KD;
vals.E3M_KD = E3M_KD;
vals.E4M_KD = E4M_KD;
vals.E3C_peak_array = E3C_peak_array;
vals.E4C_peak_array = E4C_peak_array;
vals.E3K_peak_array = E3K_peak_array;
vals.E4K_peak_array = E4K_peak_array;
vals.E3M_peak_array = E3M_peak_array;
vals.E4M_peak_array = E4M_peak_array;


if save_struct == true
save('compiled.mat', 'compiled', '-v7.3')
save('vals.mat', 'vals', '-v7.3')
else
end

% figures
if save_fig == true
% save Ket treatment figs
savefig(NaA_ket_fig, 'NaA_ket_fig')
savefig(K_ket_fig,'K_ket_fig')
savefig(NaI_ket_fig, 'NaI_ket_fig')
% save Mem treatment figs
savefig(NaA_mem_fig, 'NaA_mem_fig')
savefig(K_mem_fig,'K_mem_fig')
savefig(NaI_mem_fig, 'NaI_mem_fig')
% save Na Overlay Figs
savefig(NaC_overlay,'NaC_overlay')
savefig(NaK_overlay,'NaC_overlay')
savefig(NaM_overlay,'NaC_overlay')
else 
end


% --------------
function alignedArray = peakAlign(y,time)
%
% input arguments
%   array = x * y array output from ephysIO
%   time = [a,b] of time points in data points
%
% output arguments
%   alignedArray = peak aligned array of peak within time
%
% example usage
%   S = ephysIO(clampfile)
%   x = S.array(:,1)
%   y = S.array(:,2:end)
%   time = [1000,4000]
%   alignedArray = peakAlign(y,time)
%
% Note
%   Data loss is quite aggressive, but works for this usage. Consider
%   improving in future if use is to become more ubiquitous in code.
%% CODE

% locate minimum indexes within region of interest
 [~,ind] = min(y(time(1):time(2),:));
corrInd = ind + time(1); % correct for drift associated with peak detection

% determine the peak differences relative to the smallest
[val,~] = min(corrInd); % earliest peak and wave number
indDiff = (corrInd - val)+1; % differences in each peak index

% trim the start of the arrays, and concatenate arrays
alignedArray = zeros(sum(time)*2,size(y,2));
for i = 1:size(y,2)
    alignedArray(:,i) = y(indDiff(i):(indDiff(i)+(sum(time)*2))-1,i);
end

end

%----------------

function [] = ...
    IDplots(condition1,title1,condition2,title2)
%% Consistent variables
Vm = condition1(1).Vm;
%% Condition 1 Values
% concatenate different recordings within conditon
n = size(condition1,2);
NaID1 = zeros(size(condition1(1).NaID,2),n);
KID1 = zeros(size(condition1(1).KID,2),n);
for i = 1:n
    NaID1(:,i) = condition1(i).NaID;
    KID1(:,i) = condition1(i).KID;
end

% calculate average and SEM for each row (voltage)
NaID_mean1 = mean(NaID1,2);
KID_mean1 = mean(KID1,2);

% calculate sem for each row (voltage)
NaID_sem1 = std(NaID1,0,2)/sqrt(n);
KID_sem1 = std(KID1,0,2)/sqrt(n);

%% Condition 2 Values
% concatenate different recordings within conditon
n = size(condition2,2);
NaID2 = zeros(size(condition2(i).NaID,2),n);
KID2 = zeros(size(condition2(i).KID,2),n);
for i = 1:n
    NaID2(:,i) = condition2(i).NaID;
    KID2(:,i) = condition2(i).KID;
end

% calculate average and SEM for each row (voltage)
NaID_mean2 = mean(NaID2,2);
KID_mean2 = mean(KID2,2);

% calculate sem for each row (voltage)
NaID_sem2 = std(NaID2,0,2)/sqrt(n);
KID_sem2 = std(KID2,0,2)/sqrt(n);

%% Plotting

% set up figure environment
figure;

% plot condition one variables
NaID1_name = [title1,' Na I Density (pA/pF)'];
KID1_name = [title1,' K I Density (pA/pF)'];
errorbar(Vm,NaID_mean1,NaID_sem1,'LineStyle','none','color','black','HandleVisibility','off','linewidth',1)
hold on; errorbar(Vm,KID_mean1,KID_sem1,'LineStyle','none','color','black','HandleVisibility','off','linewidth',1)
hold on; plot(Vm,NaID_mean1,'^','color','black','DisplayName',NaID1_name,'linewidth',1.5,'MarkerSize',7,'MarkerFaceColor','white')
hold on; plot(Vm,KID_mean1,'o','color','black','DisplayName',KID1_name,'linewidth',1.5,'MarkerSize',7,'MarkerFaceColor','white')

% plot condition two variables
NaID2_name = [title2,' Na I Density (pA/pF)'];
KID2_name = [title2,' K I Density (pA/pF)'];
hold on; errorbar(Vm,NaID_mean2,NaID_sem2,'LineStyle','none','color','red','HandleVisibility','off','linewidth',1)
hold on; errorbar(Vm,KID_mean2,KID_sem2,'LineStyle','none','color','red','HandleVisibility','off','linewidth',1)
plot(Vm,NaID_mean2,'^','color','red','DisplayName',NaID2_name,'linewidth',1.5,'MarkerSize',7,'MarkerFaceColor','white')
hold on; plot(Vm,KID_mean2,'o','color','red','DisplayName',KID2_name,'linewidth',1.5,'MarkerSize',7,'MarkerFaceColor','white')

% adjust axis locations and plot labels and legends
ax = gca;
ax.XAxisLocation = 'origin'; xlabel('Membrane Potential (mV)')
ax.YAxisLocation = 'origin'; ylabel('Current Density (pA/pF)')
legend('LineWidth',1)
box off; set(gcf,'color','white'); set(gca,'linewidth',2)

end