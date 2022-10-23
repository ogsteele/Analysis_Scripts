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