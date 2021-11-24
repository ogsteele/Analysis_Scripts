function [MEA] = MEA_SpikeAnalysis(filename)
%% Author: O.G. Steele,
% Date: 06.07.18,
% Intended Use: Originally written for J.Griffiths, Cardiff Univeristy, for
%               automated analysis of multiwell systems MEA spike
%               timestamps
% Input Arguements:
%   filename = .csv exported from mutliwell systems, must be timestamps
%               with all properties
%%
Well = {'A1' 'A2' 'A3' 'A4' 'A5' 'A6' 'B1' 'B2' 'B3' 'B4' 'B5' 'B6' 'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'D1' 'D2' 'D3' 'D4' 'D5' 'D6'}';

for i = 1:24
    MEA(i) = MEA_SpikeFrequency(filename,(i-1));
    MEA(i).Well = Well(i);
end
   
end