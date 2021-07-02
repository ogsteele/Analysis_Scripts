% select the recording
recording = uigetfile('*.tdms');

figure
sgtitle('Non baseline corrected')

% whole trace
subplot(2,5,[1,2,3])
S = ephysIO(recording);
time = 4e-5*[1:size(S.array(:,2),1)]';
plot(time, S.array(:,2)*0.01)
title('whole trace')

% calculate ensemble averages
S = ephysIO(fullfile(filepath,'eventer.output/All_events/event_data.phy'));

% early
subplot(2,5,6)
early = ephysIO('pearsons_early/eventer.output/ALL_events/ensemble_average.phy');
plot(early.array(:,2))
title('early (0-600s)')

% mid
subplot(2,5,7)
mid = ephysIO('pearsons_mid/eventer.output/ALL_events/ensemble_average.phy');
plot(mid.array(:,2))
title('mid (1000-1400s)')

% late
subplot(2,5,8)
late = ephysIO('pearsons_late/eventer.output/ALL_events/ensemble_average.phy');
plot(late.array(:,2))
title('late (2000-2500s)')

% overlay
subplot(2,5,[4,5,9,10])
plot(early.array(:,2))
hold on
plot(mid.array(:,2))
plot(late.array(:,2))
title('overlay')
legend('early','mid','late')
saveas(gcf,'subtraction.pdf')

%% Penn Baseline Correction

% load in baseline.txt
baseline = table2array(readtable('20210608_005_baseline.txt')); 
early_base = median(baseline(1:60,2));
mid_base = median(baseline(100:140,2));
late_base = median(baseline(200:250,2));

% 'after' is the ephysIO structure loaded from the enesemble_average.phy from
% eventer analysis of after L689560; 'A' is the 'after' trace corrected for baseline shift

% 'before' is the ephysIO structure loaded from the enesemble_average.phy from eventer 
% analysis of before L689560; 'B' is the 'before' trace corrected for baseline shift

Early_Corrected =  early_base + early.array(:,2);
Mid_Corrected = (mid_base + mid.array(:,2)) * early_base/mid_base;
Late_Corrected = (late_base + late.array(:,2)) * early_base/late_base;

t = mid.array(:,1);
ephysIO('early.phy',[t,Early_Corrected],'S','V')
ephysIO('mid.phy',[t,Mid_Corrected],'S','V')
ephysIO('late.phy',[t,Late_Corrected],'S','V')
% Original from Teams Conversation
% t = before.array(:,1);
% 
% b = median(baseline(10:80,2)) % 'b' is the average baseline before L689560
% 
% b =
% 
%      -0.065565
% 
% a = median(baseline(end-70:end,2)) % 'a' is the average baseline after L689560
% 
% a =
% 
%     -0.06912
% 
% A = (a + after.array(:,2)) * b/a  % 'after' is the ephysIO structure loaded from the enesemble_average.phy from eventer analysis of after L689560; 'A' is the 'after' trace corrected for baseline shift
% 
% B = b + before.array(:,2) % 'before' is the ephysIO structure loaded from the enesemble_average.phy from eventer analysis of before L689560; 'B' is the 'before' trace corrected for baseline shift


%% Corrected
figure
sgtitle('Baseline corrected')

% whole trace
subplot(2,5,[1,2,3])
plot(time, S.array(:,2)*0.01)
title('whole trace')

% early
subplot(2,5,6)
plot(Early_Corrected)
title('early (0-600s)')

% mid
subplot(2,5,7)
plot(Mid_Corrected)
title('mid (1000-1400s)')

% late
subplot(2,5,8)
plot(Late_Corrected)
title('late (2000-2500s)')

% overlay
subplot(2,5,[4,5,9,10])
plot(Early_Corrected)
hold on
plot(Mid_Corrected)
plot(Late_Corrected)
title('overlay')
legend('early','mid','late')
saveas(gcf,'corrected_subtraction.pdf')