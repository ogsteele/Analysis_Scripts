function MEA = MEA_BurstFrequency(filename,WellID)
%% Author: O.G. Steele,
% Date: 06.07.18,
% Intended Use: Originally written for J.Griffiths, Cardiff Univeristy, for
%               automated analysis of multiwell systems MEA burst
%               timestamps
% Input Arguements:
%   filename = .csv exported from mutliwell systems, must be burst
%               timestamps with all properties exported
%   WellID = a number from 0-23 where 0 = A1, 23 = D6 etc (accounted for in
%            the wrapper file)
%%

MEA.Well = NaN; % Create empty template structure to default to in the absence of data

M = [21 31 12 22 32 42 13 23 33 43 24 34]'; % creates a list of channel labels 

CSV = readtable(filename); % Import CSV to table from filename input variable
RowOfInterest = find(CSV.WellID == WellID); % finds row numbers of rows containing WellID == input variable WellID

if isempty(RowOfInterest) % if no activity is detected, then populate with NaN's ...
   
    for i = 1:length(M)
        ChannelCoordinates = M(i);
            MEA.Channel(i).ChannelCoordinates = ChannelCoordinates;  
            MEA.Channel(i).Burst_Count = NaN;
            MEA.Channel(i).Burst_Frequency_in_Hz = NaN;
            MEA.Channel(i).IBI = NaN;
            MEA.Channel(i).Mean_IBI = NaN;
            MEA.Channel(i).Spike_Count = NaN;
            MEA.Channel(i).Mean_Spike_Count = NaN;
            MEA.Channel(i).Burst_Duration = NaN;
            MEA.Channel(i).Mean_Burst_Duration = NaN;
            MEA.Channel(i).Burst_Spike_Frequency = NaN;
            MEA.Channel(i).Mean_Burst_Spike_Frequency = NaN;
    end
else % ... otherwise, complete analysis
    CSVOfInterest = CSV(RowOfInterest(1):RowOfInterest(end),:); % creates new table containing only the well of interest

    for i = 1:length(M)

        ChannelCoordinates = M(i);
            MEA.Channel(i).ChannelCoordinates = ChannelCoordinates;
        ChannelNums = find(CSVOfInterest.ChannelLabel == (M(i))); % finds row numbers of rows containing ChannelLabel == 21

        if isempty(ChannelNums)
            MEA.Channel(i).Burst_Count = NaN;
            MEA.Channel(i).Burst_Frequency_in_Hz = NaN;
            MEA.Channel(i).IBI = NaN;
            MEA.Channel(i).Mean_IBI = NaN;
            MEA.Channel(i).Spike_Count = NaN;
            MEA.Channel(i).Mean_Spike_Count = NaN;
            MEA.Channel(i).Burst_Duration = NaN;
            MEA.Channel(i).Mean_Burst_Duration = NaN;
            MEA.Channel(i).Burst_Spike_Frequency = NaN;
            MEA.Channel(i).Mean_Burst_Spike_Frequency = NaN;

        else   
            ChannelofInterest = CSVOfInterest(ChannelNums(1):ChannelNums(end),:); % creates new table containing only the electrode of interest 
                                                                    % in the well of interest
            Burst_Count = height(ChannelofInterest); % returns the number of events in the channel
                MEA.Channel(i).Burst_Count = Burst_Count;
            Burst_Frequency = Burst_Count/120; % returns the frequency of events in the channel over 2 minutes
                MEA.Channel(i).Burst_Frequency_in_Hz = Burst_Frequency;
            IBI = diff(ChannelofInterest.StartTimestamp___s_); % returns the inter burst interval
                MEA.Channel(i).IBI = IBI;
            Mean_IBI = mean(diff(ChannelofInterest.StartTimestamp___s_)); % mean IBI
                MEA.Channel(i).Mean_IBI = Mean_IBI;
            Spike_Count = (ChannelofInterest.SpikeCount); % returns the spike counts from each burst
                MEA.Channel(i).Spike_Count = Spike_Count;
            Mean_Spike_Count = mean(ChannelofInterest.SpikeCount); % returns the mean spike count of all the burst in the channel
                MEA.Channel(i).Mean_Spike_Count = Mean_Spike_Count;
            Burst_Duration = (ChannelofInterest.Duration___s_); % returns the burst durations
                MEA.Channel(i).Burst_Duration = Burst_Duration;
            Mean_Burst_Duration = mean(ChannelofInterest.Duration___s_); % returns the mean burst duration from each channel
                MEA.Channel(i).Mean_Burst_Duration = Mean_Burst_Duration;
            Burst_Spike_Frequency = (ChannelofInterest.SpikeFrequency_Hz_); % returns the burst spike frequency
                MEA.Channel(i).Burst_Spike_Frequency = Burst_Spike_Frequency;
            Mean_Burst_Spike_Frequency = mean(ChannelofInterest.SpikeFrequency_Hz_); % returns the mean burst spike frequency
                MEA.Channel(i).Mean_Burst_Spike_Frequency = Mean_Burst_Spike_Frequency;
        end
    end
end
end
