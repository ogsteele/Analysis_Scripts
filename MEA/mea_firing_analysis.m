function [channel,crossCorrMatNorm,syncIndexMat,tau,syncIndex] = mea_firing_analysis(csvfilename)
% mea_firing_analysis.m
%
% Calculates some properties of cell firing extracted from MCS multi-well
% MEA system. Designed initially for Christian Schnell/Paul Kemp.
%
% INPUTS:
%   csvfilename = char array containg the filename of the .csv file
%                 exported from MCS MEA system
%
% OUTPUTS:
%   channel = struc containing the channel-by-channel output from the
%               script.
%   crossCorrMatNorm = cell array containing normalised pair-wise
%                       crosscorrelation matrices
%   syncIndexMat = matrix containing the average peak normalised
%                       crosscorrelation scores
%   tau = time lag axis for individual cross correlations
%   syncIndex = the mean synchronisation index for the well
%
% by Jon Brown, University of Exeter, 10-3-2017

%% define some analysis parameters

%%% minimum number of spikes to do analysis
minNumSpikes = 30;

%%% bins for histogram
nBins = 200;                                            % number of bins
maxLogISI = 3;                                          % max log ISI in log(s)
minLogISI = -3;                                         % min log ISI in log(s)
histBinEdges = logspace(minLogISI,maxLogISI,nBins+1);   % create bin edges

%%% complex spike index ISI threshold
complexSpikeThresh = 0.02;                             % threshold in s

%%% cross corr
ccorr_win       = [-2,2];
ccorr_method    = 1;
ccorr_NBins     = 100;
ccorr_bin_range = [-0.2,0.2];                           % ccorr range for average

%%% set up line format for raster
LineFormat.Color = [0.6,0.6,0.6];
LineFormat.LineWidth = 0.2;
LineFormat.LineStyle = '-';
VertSpikeHeight = 0.8;


%% load in .csv file
[ChannelID,...
    ChannelLabel,...
    WellID,...
    WellLabel,...
    CompoundID,...
    CompoundName,...
    Experiment,...
    DosepM,...
    DoseLabel,...
    Timestamps,...
    MaximumAmplitudepV,...
    MinimumAmplitudepV,...
    PeaktopeakAmplitudepV] = import_csv_mea(csvfilename);

%% convert Timestamps to seconds - currently microsec
Timestamps = Timestamps*1e-6;

%% organise into different channels

[C,ia,ic] = unique(ChannelID,'stable');

% number of unique channels
numCh = length(C);

% make a struct of each individual channel

% preallocate channel struc
channel(numCh).ChannelID = [];
channel(numCh).ChannelLabel = [];
channel(numCh).WellID = [];
channel(numCh).WellLabel = [];
channel(numCh).CompoundID = [];
channel(numCh).CompoundName = [];
channel(numCh).Experiment = [];
channel(numCh).numSpikes = [];
channel(numCh).PeaktopeakAmplitudepV = [];
channel(numCh).Timestamps = [];
channel(numCh).ISI = [];
channel(numCh).ISIhist = [];
channel(numCh).complexSpikeIndex = [];

for i = 1:numCh
    % add metadata to struc
    channel(i).ChannelID        = ChannelID(ia(i));
    channel(i).ChannelLabel     = ChannelLabel(ia(i));
    channel(i).WellID           = WellID(ia(i));
    channel(i).WellLabel        = WellLabel(ia(i));
    channel(i).CompoundID       = CompoundID(ia(i));
    channel(i).CompoundName     = CompoundName(ia(i));
    channel(i).Experiment       = Experiment(ia(i));
    
    % add timestamps and amplitude to struc
    if i<numCh
        channel(i).Timestamps               = Timestamps(ia(i):ia(i+1)-1);
        channel(i).PeaktopeakAmplitudepV    = PeaktopeakAmplitudepV(ia(i):ia(i+1)-1);
    else
        channel(i).Timestamps               = Timestamps(ia(i):end);
        channel(i).PeaktopeakAmplitudepV    = PeaktopeakAmplitudepV(ia(i):end);
    end
    channel(i).numSpikes = length(channel(i).Timestamps);
end

%% plot raster plots
figure
hold on
% put timestamps in correct format for plotSpikeRater.m
spikesT = cell(numCh,1);
for i = 1:numCh
    spikesT{i} = channel(i).Timestamps';
end
plotSpikeRaster(spikesT,'PlotType','vertline',...
    'LineFormat',LineFormat,...
    'VertSpikeHeight',VertSpikeHeight);
xlabel('time (s)')
drawnow

%% calculate ISI histogram for each cell

for i = 1:numCh
    % only caclulate complex spike index if channel has a threshold number
    % of spikes, defined by minNumSpikes
    if channel(i).numSpikes>=minNumSpikes
        channel(i).ISI = diff(channel(i).Timestamps);                               % calculate ISI (in s)
        [channel(i).ISIhist.N,channel(i).ISIhist.binEdges] = ...
            histcounts(channel(i).ISI,histBinEdges,'Normalization','probability');  % calculate ISI hist
        channel(i).ISIhist.binCentres = diff(channel(i).ISIhist.binEdges) + channel(i).ISIhist.binEdges(1:end-1);
        
        ind = find(channel(i).ISIhist.binEdges <= complexSpikeThresh,1,'last');     % find bin edge
        channel(i).complexSpikeIndex = sum(channel(i).ISIhist.N(1:ind-1));          % proportion of spikes below ISI threshold
        
             figure
             semilogx(channel(i).ISIhist.binCentres,channel(i).ISIhist.N)
    end
end

%% calculate cross corr matrix

% only for channels with a a threshold number of spikes, defined by 
% minNumSpikes
numSpikesPerCh = zeros(numCh,1);
for i=1:numCh
    numSpikesPerCh(i) = channel(i).numSpikes;
end
indCh = find(numSpikesPerCh>=minNumSpikes);
numChRed = length(indCh);

% preallocate
crossCorrMat        = cell(numChRed);
crossCorrMatNorm    = cell(numChRed);
syncIndexMat = zeros(numChRed);

for i = 1:numChRed
    for j = 1:numChRed
        [crossCorrMat{i,j},tau] = ccorr(channel(indCh(i)).Timestamps,...
            channel(indCh(j)).Timestamps,...
            ccorr_win,...
            'n',...
            [],...
            ccorr_method,...
            0,...
            ccorr_NBins);
        % normalise the crosscorrelation
        % 1) divide by the sum of all the bins
        % 2) multiply by the number of bins
        % 3) subtract 1
        crossCorrMatNorm{i,j} = ((crossCorrMat{i,j}/sum(crossCorrMat{i,j}))*ccorr_NBins)-1;
        % find the tau range of interest to calculate synch index
        ind = tau>=ccorr_bin_range(1) & tau<=ccorr_bin_range(2);
        syncIndexMat(i,j) = mean(crossCorrMatNorm{i,j}(ind));
    end
end

% remove the indentity matrix values (i.e. autocorrelations)
I = logical(eye(numChRed));
syncIndexMat(I) = NaN;

%% calculate the mean synchronization index; ~0 is no synchrony >0 more synchronised
% extract upper triangular part of syncIndexMat
ind = logical(triu(ones(numChRed),1));
% mean of the upper  triangle - the lower triangle is mirror image
syncIndex = mean(syncIndexMat(ind));

%% plot the cross corr matrix
figure
imagesc(syncIndexMat)
set(gca,'YDir','normal','Box','off')
xlabel('Channel #')
ylabel('Channel #')
colormap(jet)
caxis([0,4])
title(['mean synchronicity index = ' num2str(syncIndex,2)])
colorbar

end
