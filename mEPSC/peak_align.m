E3M_NMDAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E3M/NMDAR/E3M_NMDAR_average.h5');
E3F_NMDAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E3F/NMDAR/E3F_NMDAR_average.h5');
E4M_NMDAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E4M/NMDAR/E4M_NMDAR_average.h5');
E4F_NMDAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E4F/NMDAR/E4F_NMDAR_average.h5');


%% Load and concatenate

% add time as the first and save as ensemble average
NMDAR_aves = [E3M_NMDAR.array(:,2) E3F_NMDAR.array(:,2) E4M_NMDAR.array(:,2) E4F_NMDAR.array(:,2)];
time = 5e-5*[1:size(NMDAR_aves,1)]';
% create the output in the structure
filename = 'NMDAR_ensemble_averages.phy';
% save the ensemble average
ephysIO(char(filename),[time NMDAR_aves],'s','A') 

%%
S = ephysIO('NMDAR_peak_scaled.h5');
array = S.array(:,2:end);
% apply a filter to clean the look of the data
YF = filter1(array, time, 0, 500);
figure;
plot(YF,'linewidth',2)
legend('E3M','E3F','E4M','E4F')


figure
plot(S.array(:,2),'linewidth',0.5,'color',[0.9290,0.6940, 0.1250 0.1])
hold on
plot(smooth((median(S.array(:,2:end-1),2))),'linewidth',2,'color',[0.9290,0.6940, 0.1250 1])
box off
set(gcf,'color','w');
set(gca,'visible','off')



%% 
E3M_AMPAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E3M/AMPAR/E3M_AMPAR_average.h5');
E3F_AMPAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E3F/AMPAR/E3F_AMPAR_average.h5');
E4M_AMPAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E4M/AMPAR/E4M_AMPAR_average.h5');
E4F_AMPAR = ephysIO('/Volumes/T7/Experiments/NMDAR mEPSC/04_processing/Subtraction/Averages/E4F/AMPAR/E4F_AMPAR_average.h5');

% add time as the first and save as ensemble average
AMPAR_aves = [E3M_AMPAR.array(:,2) E3F_AMPAR.array(:,2) E4M_AMPAR.array(:,2) E4F_AMPAR.array(:,2)];
time = 5e-5*[1:size(AMPAR_aves,1)]';
% create the output in the structure
filename = 'AMPAR_ensemble_averages.phy';
% save the ensemble average
ephysIO(char(filename),[time AMPAR_aves],'s','A') 

S = ephysIO('AMPAR_peak_scaled.h5');
array = S.array(:,2:end);
% apply a filter to clean the look of the data
YF = filter1(array, time, 0, 500);
figure;
plot(YF,'linewidth',2)
set(gcf,'color','w');
set(gca,'visible','off')
legend('E3M','E3F','E4M','E4F')

figure
plot(S.array(:,2),'linewidth',0.5,'color',[0.9290,0.6940, 0.1250 0.1])
hold on
plot(smooth((median(S.array(:,2:end-1),2))),'linewidth',2,'color',[0.9290,0.6940, 0.1250 1])
box off
set(gcf,'color','w');
set(gca,'visible','off')
