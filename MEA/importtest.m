function Out = importtest(csvfilename)
%% Options for Analysis

% 1 = Yes
% 0 = No

CrossCorrPlot = 0;
RasterPlot = 0;
ThreeDimensional = 0;

%% Filename operator

Out.filename = csvfilename;
%% Operational Parameters

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


%% Locate index of well positions
[a,b,c] = import_csv_mea(csvfilename);

[C,ia,ic] = unique(c,'stable');

for x = 1:((length(C))-1)
        ia(x,2) = ia((x+1),1) - 1;
end
    ia(length(C),2) = length(ic);
    


%% load in .csv file of interest (seperated into wells)

for x = 1:length(C)
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
    PeaktopeakAmplitudepV] = import_csv_mea(csvfilename,((ia(x,1))+1),((ia(x,2))+1));

%% convert Timestamps to seconds - currently microsec
Timestamps = Timestamps*1e-6;

%% organise into different channels

[C2,ia2,ic2] = unique(ChannelID,'stable');

% number of unique channels
numCh = length(C2);

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
    Out.Well(x).channel(i).ChannelID        = ChannelID(ia2(i));
    Out.Well(x).channel(i).ChannelLabel     = ChannelLabel(ia2(i));
    Out.Well(x).channel(i).WellID           = WellID(ia2(i));
    Out.Well(x).channel(i).WellLabel        = WellLabel(ia2(i));
    Out.Well(x).channel(i).CompoundID       = CompoundID(ia2(i));
    Out.Well(x).channel(i).CompoundName     = CompoundName(ia2(i));
    Out.Well(x).channel(i).Experiment       = Experiment(ia2(i));
    
    % add timestamps and amplitude to struc
    if i<numCh
        Out.Well(x).channel(i).Timestamps               = Timestamps(ia2(i):ia2(i+1)-1);
        Out.Well(x).channel(i).PeaktopeakAmplitudepV    = PeaktopeakAmplitudepV(ia2(i):ia2(i+1)-1);
    else
        Out.Well(x).channel(i).Timestamps               = Timestamps(ia2(i):end);
        Out.Well(x).channel(i).PeaktopeakAmplitudepV    = PeaktopeakAmplitudepV(ia2(i):end);
    end
    Out.Well(x).channel(i).numSpikes = length(Out.Well(x).channel(i).Timestamps);
end
%% calculate ISI histogram for each cell

for i = 1:numCh
    % only caclulate complex spike index if channel has a threshold number
    % of spikes, defined by minNumSpikes
    if Out.Well(x).channel(i).numSpikes>=minNumSpikes
        Out.Well(x).channel(i).ISI = diff(Out.Well(x).channel(i).Timestamps);                               % calculate ISI (in s)
        [Out.Well(x).channel(i).ISIhist.N,Out.Well(x).channel(i).ISIhist.binEdges] = ...
            histcounts(Out.Well(x).channel(i).ISI,histBinEdges,'Normalization','probability');  % calculate ISI hist
        Out.Well(x).channel(i).ISIhist.binCentres = diff(Out.Well(x).channel(i).ISIhist.binEdges) + Out.Well(x).channel(i).ISIhist.binEdges(1:end-1);
        
        ind = find(Out.Well(x).channel(i).ISIhist.binEdges <= complexSpikeThresh,1,'last');     % find bin edge
        Out.Well(x).channel(i).complexSpikeIndex = sum(Out.Well(x).channel(i).ISIhist.N(1:ind-1));          % proportion of spikes below ISI threshold
        
             %figure
             %semilogx(channel(i).ISIhist.binCentres,channel(i).ISIhist.N)
    end
end

%% plot raster plots

if RasterPlot == 1
    
    figure
    set(gcf,'color','w');
    subplot(2,2,[1,2]);
    hold on
    % put timestamps in correct format for plotSpikeRater.m
    spikesT = cell(numCh,1);
    for i = 1:numCh
        spikesT{i} = Out.Well(x).channel(i).Timestamps';
    end
    plotSpikeRaster(spikesT,'PlotType','vertline',...
        'LineFormat',LineFormat,...
        'VertSpikeHeight',VertSpikeHeight);
    xlabel('time (s)')
    ylabel('Electrode number')
    drawnow
end

%% calculate cross corr matrix

% only for channels with a a threshold number of spikes, defined by 
% minNumSpikes
numSpikesPerCh = zeros(numCh,1);
for i=1:numCh
    numSpikesPerCh(i) = Out.Well(x).channel(i).numSpikes;
end
indCh = find(numSpikesPerCh>=minNumSpikes);
numChRed = length(indCh);

% preallocate
crossCorrMat        = cell(numChRed);
crossCorrMatNorm    = cell(numChRed);
syncIndexMat = zeros(numChRed);

for i = 1:numChRed
    for j = 1:numChRed
        [crossCorrMat{i,j},tau] = ccorr(Out.Well(x).channel(indCh(i)).Timestamps,...
            Out.Well(x).channel(indCh(j)).Timestamps,...
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
Out.syncIndex(x) = nanmean(syncIndexMat(ind));

%% plot the cross corr matrix
if CrossCorrPlot == 1
    subplot(2,2,3);
    imagesc(syncIndexMat)
    set(gca,'YDir','normal','Box','off')
    xlabel('Channel #')
    ylabel('Channel #')
    colormap(jet)
    caxis([0,4])
    title(['mean synchronicity index = ' num2str(syncIndex,2)])
    colorbar
end

%% 3D Bar Plots
   
if ThreeDimensional == 1

    % creates array where collumn 1 is el pos, 2 is spike num
    SpikeNums = [Out.Well(x).channel(:).numSpikes];
    ElPos = [Out.Well(x).channel(:).ChannelLabel];
    ElPos_SpikeNum = [(ElPos)',(SpikeNums)'];
    
    Corrected_ElPos_SpikeNum = zeros(12,2);
    Corrected_ElPos_SpikeNum(:,1) = [21,31,12,22,32,42,13,23,33,43,24,34]';
    
    El21 = find(ElPos_SpikeNum(:,1) == 21);
        if El21 == []
            Corrected_ElPos_SpikeNum(1,2) = 0;
        else
            Corrected_ElPos_SpikeNum(1,2) = ElPos_SpikeNum(El21,2);
        end
    El31 = find(ElPos_SpikeNum(:,1) == 31);
        if El31 == []
            Corrected_ElPos_SpikeNum(2,2) = 0;
        else
            Corrected_ElPos_SpikeNum(2,2) = ElPos_SpikeNum(El31,2);
        end
    EL12 = find(ElPos_SpikeNum(:,1) == 12);
        if El12 == []
            Corrected_ElPos_SpikeNum(3,2) = 0;
        else
            Corrected_ElPos_SpikeNum(3,2) = ElPos_SpikeNum(El12,2);
        end
    EL22 = find(ElPos_SpikeNum(:,1) == 22);
        if El22 == []
            Corrected_ElPos_SpikeNum(4,2) = 0;
        else
            Corrected_ElPos_SpikeNum(4,2) = ElPos_SpikeNum(El22,2);
        end
    EL32 = find(ElPos_SpikeNum(:,1) == 32);
        if El32 == []
            Corrected_ElPos_SpikeNum(5,2) = 0;
        else
            Corrected_ElPos_SpikeNum(5,2) = ElPos_SpikeNum(El32,2);
        end
    El42 = find(ElPos_SpikeNum(:,1) == 42);
        if El42 == []
            Corrected_ElPos_SpikeNum(6,2) = 0;
        else
            Corrected_ElPos_SpikeNum(6,2) = ElPos_SpikeNum(El42,2);
        end
    El13 = find(ElPos_SpikeNum(:,1) == 13);
        if El13 == []
            Corrected_ElPos_SpikeNum(7,2) = 0;
        else
            Corrected_ElPos_SpikeNum(7,2) = ElPos_SpikeNum(El13,2);
        end
    El23 = find(ElPos_SpikeNum(:,1) == 23);
        if El23 == []
            Corrected_ElPos_SpikeNum(8,2) = 0;
        else
            Corrected_ElPos_SpikeNum(8,2) = ElPos_SpikeNum(El23,2);
        end
    El33 = find(ElPos_SpikeNum(:,1) == 33);
        if El33 == []
            Corrected_ElPos_SpikeNum(9,2) = 0;
        else
            Corrected_ElPos_SpikeNum(9,2) = ElPos_SpikeNum(El33,2);
        end
    El43 = find(ElPos_SpikeNum(:,1) == 43);
        if El43 == []
            Corrected_ElPos_SpikeNum(10,2) = 0;
        else
            Corrected_ElPos_SpikeNum(10,2) = ElPos_SpikeNum(El43,2);
        end
    El24 = find(ElPos_SpikeNum(:,1) == 24);
        if El24 == []
            Corrected_ElPos_SpikeNum(11,2) = 0;
        else
            Corrected_ElPos_SpikeNum(11,2) = ElPos_SpikeNum(El24,2);
        end
    El34 = find(ElPos_SpikeNum(:,1) == 34);
        if El34 == []
            Corrected_ElPos_SpikeNum(12,2) = 0;
        else
            Corrected_ElPos_SpikeNum(12,2) = ElPos_SpikeNum(El34,2);
        end
    
    
    spikecounts = zeros(4,4)
    spikecounts = [0,
    (Corrected_ElPos_SpikeNum(1,2)),
    (Corrected_ElPos_SpikeNum(2,2)),
    0;...
    (Corrected_ElPos_SpikeNum(3,2)),
    (Corrected_ElPos_SpikeNum(4,2)),
    (Corrected_ElPos_SpikeNum(5,2)),
    (Corrected_ElPos_SpikeNum(6,2));...
    (Corrected_ElPos_SpikeNum(7,2)),
    (Corrected_ElPos_SpikeNum(8,2)),
    (Corrected_ElPos_SpikeNum(9,2)),
    (Corrected_ElPos_SpikeNum(10,2));...
    0,
    (Corrected_ElPos_SpikeNum(11,2)),
    (Corrected_ElPos_SpikeNum(12,2)),
    0];


%spikecounts =  [0,(Well(x).channel(1).numSpikes),(Well(x).channel(2).numSpikes),0;...
                %(Well(x).channel(3).numSpikes),(Well(x).channel(4).numSpikes),(Well(x).channel(5).numSpikes),(Well(x).channel(6).numSpikes);...
                %(Well(x).channel(7).numSpikes),(Well(x).channel(8).numSpikes),(Well(x).channel(9).numSpikes),(Well(x).channel(10).numSpikes);...
                %0,(Well(x).channel(11).numSpikes),(Well(x).channel(12).numSpikes),0];

subplot(2,2,4);
b = bar3(spikecounts);
colorbar

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

ylabel('Y Electrode Coordinate')
xlabel('X Electrode coordinate')
zlabel('Spike Count')
else 
end
end
end