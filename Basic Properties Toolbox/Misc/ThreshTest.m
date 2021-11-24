function [Thresh] = ThreshTest(filename1,sweep)
CS = abfload(filename1);
CSs = CS(:,sweep);
[peakmv,peaktime] = max(CSs);
    PeakmV = peakmv

tdiffCSs = smooth(diff(diff(diff(CSs))));
[ThreshPeak,ThreshInd] = max(tdiffCSs);
Thresh = CSs(ThreshInd);


%%  Find min of diff(diff(diff(CSs), then find index of local max and apply index to raw trace with AP to find threshold