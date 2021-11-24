function [PeakmV,AftHyp,SpikeAmp] = SpikeAnalysis(SpikeSweepWindow,pre,peak,post)
PeakmV = max(SpikeSweepWindow(pre:post,:,:));
    %xlswrite(SpreadSheet,PeakAmp,'Sheet1','B31');
AftHyp = min(SpikeSweepWindow(peak:post,:,:));
SpikeAmp = PeakmV - AftHyp;
end 