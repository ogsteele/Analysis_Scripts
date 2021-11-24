function MEA = MEA_SpikeFrequency(filename,WellID)
%% Author: O.G. Steele,
% Date: 06.07.18,
% Intended Use: Originally written for J.Griffiths, Cardiff Univeristy, for
%               automated analysis of multiwell systems MEA spike
%               timestamps
% Input Arguements:
%   filename = .csv exported from mutliwell systems, must be spike
%               timestamps with all properties exported
%   WellID = a number from 0-23 where 0 = A1, 23 = D6 etc (accounted for in
%            the wrapper file
%%
MEA.Well = NaN; % Create empty template structure to default to in the absence of data

M = [21 31 12 22 32 42 13 23 33 43 24 34]'; % creates a list of channel labels 

CSV = readtable(filename); % Import CSV to table from filename input variable
RowOfInterest = find(CSV.WellID == WellID); % finds row numbers of rows containing WellID == input variable

if isempty(RowOfInterest) % if no activity is detected, then populate with NaN's ...
   
    for i = 1:length(M)
        ChannelCoordinates = M(i);
            MEA.Channel(i).ChannelCoordinates = ChannelCoordinates;  
            MEA.Channel(i).SpikeCount = NaN;
            MEA.Channel(i).SpikeFrequency = NaN;
            MEA.Channel(i).ISI = NaN;
            MEA.Channel(i).MeanISI = NaN;
            MEA.Channel(i).MeanP2PAmplitude = NaN;
    end
else       % ... otherwise, complete analysis
    CSVOfInterest = CSV(RowOfInterest(1):RowOfInterest(end),:); % creates new table containing only the well of interest

    for i = 1:length(M)

        ChannelCoordinates = M(i);
            MEA.Channel(i).ChannelCoordinates = ChannelCoordinates;
        ChannelNums = find(CSVOfInterest.ChannelLabel == (M(i))); % finds row numbers of rows containing ChannelLabel == 21

        if isempty(ChannelNums) % if no activity is detected in the electrodes, populate with NaN's
            MEA.Channel(i).SpikeCount = NaN;
            MEA.Channel(i).SpikeFrequency = NaN;
            MEA.Channel(i).ISI = NaN;
            MEA.Channel(i).MeanISI = NaN;
            MEA.Channel(i).MeanP2PAmplitude = NaN;

        else   
            ChannelofInterest = CSVOfInterest(ChannelNums(1):ChannelNums(end),:); % creates new table containing only the electrode of interest 
                                                                    % in the well of interest
            SpikeCount = height(ChannelofInterest); % returns the number of events in the channel
                MEA.Channel(i).SpikeCount = SpikeCount;
            SpikeFrequency = SpikeCount/120; % returns the frequency of events in the channel
                MEA.Channel(i).SpikeFrequency = SpikeFrequency;
            ISI = diff(ChannelofInterest.Timestamp___s_); % returns the inter spike interval
                MEA.Channel(i).ISI = ISI;
            MeanISI = mean(diff(ChannelofInterest.Timestamp___s_)); % mean ISI
                MEA.Channel(i).MeanISI = MeanISI;
            MeanP2PAmplitude = mean(ChannelofInterest.Peak_to_peakAmplitude_pV_); % returns the mean of the P2P amplitudes
                MEA.Channel(i).MeanP2PAmplitude = MeanP2PAmplitude;
        end
    end
end

