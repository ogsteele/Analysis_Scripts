function scaled_traces = peakscale(array)
%% peakscale
% scale traces to their peak

% input should be an array of traces with the first collumn as time in seconds, eg out
% from S = ephysIO

% following would save the corrected together
%   Comp = ephysIO('corr_ens_ave_compound.phy')
%   Ampa = ephysIO('corr_ens_ave_AMPA.phy')
%   Nmda = ephysIO('corr_ens_ave_NMDA.phy')
%   array = [Comp.array Ampa.array(:,2) Nmda.array(:,2)]
%   ephysIO('conc_test.phy',array,'S','V')

%% Example
% read in ensemble averages with ephysIO if necessary then run function
% S = ephysIO('conc_test.phy');
% array = S.array;
% scaled_traces = peakscale(array);
%% plot raw
%% plot raw and save too
figure; plot(array(:,1),array(:,2:end));
title('Raw Traces'); xlabel('Time (s)'); ylabel('Voltage (V)')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
legend('Compound','AMPAR','NMDAR')
saveas(gcf,'raw.pdf')
saveas(gcf,'raw.emf')

%% measure baseline
% plot traces overlayed
figure; plot(array(:,2:end)); 
title('1. Select baseline region'); xlabel('Data points'); ylabel('Voltage (V)')
% select and plot baseline region
[b,~] = ginput(2);
hold on; xline(b(1),'color','green','linestyle','--')
hold on; xline(b(2),'color','green','linestyle','--')
% get mean baseline values
round_b = floor(b);
mean_baseline = zeros(1,(size(array,2)-1));
for i = 2:size(array,2)
    mean_baseline(i-1) = mean(array(round_b(1):round_b(2),i));
end
%% subtract baselines
for i = 2:size(array,2)
    array(:,i) = (array(:,i) - mean_baseline(i-1));
end
%% measure peak amplitudes
% select the baseline region
title('2. Select peak region')
[p,~] = ginput(2);
hold on; xline(p(1),'color','red','linestyle','--')
hold on; xline(p(2),'color','red','linestyle','--')
% get mean peak values (from baseline subtracted trace)
round_p = floor(p);
max_peak = zeros(1,(size(array,2)-1));
mean_peak = zeros(1,(size(array,2)-1));
max_ind = zeros(1,(size(array,2)-1));
for i = 2:size(array,2)
    [max_peak(i-1),max_ind(i-1)] = max(array(round_p(1):round_p(2),i));
end
max_ind = max_ind + round_p(1);
for i = 2:size(array,2)
    mean_peak(i-1) = mean(array(max_ind(i-1)-10:max_ind(i-1)+10,i));
end
% plot peak values (back onto non baseline subtracted trace)
for i = 1:size(mean_peak,2)
    hold on; plot(max_ind(i),mean_peak(i)+mean_baseline(i),'ro','MarkerSize',10,'linewidth',3)
end

%% Scale the traces
% Calculate scale factor to make peak equal to the mean peak amplitude
scale_factor = mean_peak / mean(mean_peak);
% Scale the traces and apply offset equal to the mean baseline
scaled_traces = array;
for i = 2:size(array,2)
    scaled_traces(:,i) = array(:,i) / scale_factor(i-1) + mean(mean_baseline);
end
% plot peak scaled traces
figure; plot(scaled_traces(:,1),scaled_traces(:,2:end)); box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Voltage (V)'); title('Peak Scaled traces')
legend('Compound','AMPAR','NMDAR')
saveas(gcf,'peakscaled.pdf')
saveas(gcf,'peakscaled.emf')
%% Save outputs
% figure
%saveas(gcf,'peakscaled.mfig')
% data
ephysIO('peakscaled.phy',scaled_traces,'S','A')
