% select cell of interest
close all; clear

path = uigetdir(); cd(path);
myFolders = split(path,"\");
%details = char(fullfile(myFolders(6),myFolders(7),...
%                        myFolders(8),myFolders(9),...
%                        myFolders(10),myFolders(11)));
%subt = ['...\',details];
clipboard('copy',path) % copy to the clipboard

% load raw events and plot average
M=ephysIO('mlm_before/eventer.output/ALL_events/event_data.phy');
time = M.array(:,1);
m_event = M.array(:,2:end)*1e12;
%figure; set(gcf, 'Position',  [100, 100, 1200, 400]); subplot(1,3,1)
figure;
plot(time,m_event,'color',[0.5 0.5 0.5 0.1],'HandleVisibility','off')
ylabel('Amplitude (pA)'); xlabel('Time (s)');
hold on; plot(time,median(m_event,2), 'color', 'blue', 'linewidth' ,2)
m_ensemble = median(m_event,2);
legend('Ensemble Average','linewidth',1)
xlim([0 0.1]);ylim([-60 30])
box off; set(gca,'linewidth',2); set(gcf,'color','white'); title('Mixed')

% load AMPAR events and plot average
A=ephysIO('mlm_after/eventer.output/ALL_events/event_data.phy');
a_event = A.array(:,2:end)*1e12;
%subplot(1,3,2); 
figure;
plot(time,a_event,'color',[0.5 0.5 0.5 0.1],'HandleVisibility','off')
ylabel('Amplitude (pA)'); xlabel('Time (s)');
hold on; plot(time,median(a_event,2), 'color', 'red', 'linewidth' ,2)
a_ensemble = median(a_event,2);
legend('Ensemble Average','linewidth',1)
xlim([0 0.1]);ylim([-60 30])
box off; set(gca,'linewidth',2); set(gcf,'color','white'); title('AMPAR')

% subtract the values
%subplot(1,3,3); 
figure; 
plot(time,m_ensemble, 'color', 'blue', 'linewidth' ,2)
hold on; plot(time,a_ensemble, 'color', 'red', 'linewidth' ,2)
n_ensemble = m_ensemble-a_ensemble;
hold on; plot(time,n_ensemble, 'color', 'yellow', 'linewidth' ,2);
legend('Mixed Ensemble', 'AMPAR Ensemble', 'NMDAR Ensemble', 'linewidth',1);
xlim([0 0.1]);ylim([-40 20])
box off; set(gca,'linewidth',2); set(gcf,'color','white'); title('Subtraction')
ylabel('Amplitude (pA)'); xlabel('Time (s)');
warning('off')
%sgtitle(subt);

% saving
mkdir output_2022; cd('output_2022') % make new directory and change to it
saveas(gcf,'subtraction_2022.pdf');saveas(gcf,'subtraction_2022') % save subtraction figure
ephysIO('AMPAR_ensemble.phy',[time a_ensemble*1e-12],'s','A')
ephysIO('AMPAR_events.phy',[time a_event*1e-12],'s','A')
ephysIO('mixed_ensemble.phy',[time m_ensemble*1e-12],'s','A')
ephysIO('mixed_events.phy',[time m_event*1e-12],'s','A')
ephysIO('NMDAR_ensemble.phy',[time n_ensemble*1e-12],'s','A')


