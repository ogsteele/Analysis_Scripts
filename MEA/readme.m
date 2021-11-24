%% MEA_Spike/BurstAnalysis Documentation
% Author: O.G. Steele
% Date: 06.07.18

%% Calling Wrappers
% The wrappers, MEA_BurstAnalysis and MEA_SpikeAnalysis, call the functions
% MEA_BurstFrequency and MEA_SpikeFrequency respectively repeating for the
% number of wells in the plate (24 in this instance.)
% To produce the structure x containing the rearranged and calculated .csv timestamp
% file a, input the following command:
% x = MEA_BurstAnalysis('a.csv'); or x = MEA_SpikeAnalysis('a.csv');
% Please ensure all folders are in the correct path

%% Splitting structures
% To split the structure x into a smaller structure a containing the
% 1st, 2nd, 3rd and nth element of structure x, input the following 
% commannd:
% a = [x(1),x(2),x(3),x(n)];
% For example, KO = [Spike_MEA(1),Spike_MEA(4),Spike_MEA(6:9)] would
% contain the 1st, 4th, 6th, 7th, 8th and 9th rows of the structure 'Spike_MEA'