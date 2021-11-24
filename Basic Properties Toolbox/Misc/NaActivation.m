function [NaAPeaks,NaADensity,NaADrivingForce,NaAConductance] = NaActivation (filename3,capacitance)
% CONSTANTS
pA = linspace(-120,60,37)';
DF = pA-66.68;
% Prerequisites
FileTemp = abfload(filename3);
DavidActInactPro = reshape(FileTemp,10322,37);
% Calculations
NaAVal = DavidActInactPro(2924:5001,:,:);
NaAPeaks = (min(NaAVal))';
NaADensity = NaAPeaks/capacitance;
NaADrivingForce = NaADensity./DF;
NaAConductance = pA./max(NaADrivingForce);
end
