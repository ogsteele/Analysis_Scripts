function output = ACQ4_IO_Batch
% Date of Creation: 21.02.19
% Updated: 21.02.19

%% Parameters
% Define Sample Rate in kHz
samplerate = 40;


%% Analysis
% Clear windows
close all

% Define a starting folder.
cd('/Users/ogsteele/OneDrive - University of Sussex/MATLAB/PhD');

% Get xlsx spreadsheet with filenames of top directories of I/V curves
[fileName,dirName]=uigetfile('.xlsx');
[~,~,rawData] = xlsread(fullfile(dirName,fileName));

% Define for loop iterations 
for n = 1:length(rawData)
    

    % Ask user to confirm or change.
    topLevelFolder = char(rawData(n));

    % Get list of all subfolders.
    allSubFolders = genpath(topLevelFolder);

    % Parse into a cell array.
    remain = allSubFolders;
    listOfFolderNames = {};
    while true
        [singleSubFolder, remain] = strtok(remain, ':');
        if isempty(singleSubFolder)
            break;
        end
        listOfFolderNames = [listOfFolderNames singleSubFolder];
    end
    numberOfFolders = length(listOfFolderNames);

    % Change directory to first sweep 
    % Note, not first folder as that's the master folder
    for i = 2:numberOfFolders
        cd(char(listOfFolderNames(:,i)));
        Data = h5read('Clamp2.ma','/data');
        Trace(:,i-1) = Data(:,1);
    end

    % adjust sweeps to zero and smooth to compensate for noise
    basemean = mean(Trace(1:3500,:));
    adjusted = Trace(:,:) - basemean;

    % plot the sweeps and get points pre and post field potential
    plot(adjusted(3500:5500,:));
    ylim([-0.004 0.005]);
    xlim([0 2000]);
    hold on
    [x,y] = getpts;
    close

    % Seperate out window of interest
    fEPSP_start = round(x(1)+3500);
    fEPSP_end = round(x(2)+3500);
    fEPSP_Window = adjusted((fEPSP_start:fEPSP_end),:);

    % Calculate time of window
    Time_Conv = length(fEPSP_Window)/samplerate;
    xtime = linspace(0,Time_Conv,length(fEPSP_Window));

    % Plot fEPSP Traces
    %subplot(1,2,1)
    figure
    for l = 1:10
        plot(xtime,(smooth(smooth((fEPSP_Window(:,l)*1000)))))
        hold on
    end

    ylabel('Amplitude (mV)');
    xlabel('Time (ms)');
    xlim([0 Time_Conv]);
    legend({'0 V','10 V','20 V','30 V','40 V','50 V','60 V','70 V','80 V','90 V'},'Location','southeast');
    title('fEPSP I/V Traces')
    hold off

    % Index of maximum amplitude, and incidentally the peak amplitude of each
    % stimulation intensity
    for o = 1:10
        [M(o),I(o)] = min((fEPSP_Window(:,o)));
    end
    M = M*1000;

    % index of 30% and 70% of maximum amplitudes
    I_slope_30 = round(I*0.3);
    I_slope_70 = round(I*0.7);

    % calculate fEPSP slope between 30% and 70% peak amplitude
    for k = 1:10
        slope(k) = mean(gradient(fEPSP_Window((I_slope_30(k):I_slope_70(k)),k)));
    end

    % Half Maximal Stimulation Intensity
    StimIntensity = linspace(0,90,10)';
    HM = interp1(M,StimIntensity,(min(M)/2));
    MaxAmp = min(M);

    % Plot I/O Curve
    % subplot(1,2,2)
    figure
    plot(StimIntensity,M');
    hold on
    ylabel('Amplitude (mV)');
    xlabel('Stimulation Intensity (V)');
    title('I/O Curve');

    % Annotate with Half Maximal Stimulation Intensity
    text = 'Half maximal stimulation intensity (V) = %f';
    str = sprintf(text,HM);
    dim = [0.441071428571428 0.832142857142859 0.43125 0.053571428571429];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');

    % Create output structure
    output.MaxAmp(n) = MaxAmp;
    output.HalfMax(n) = HM;
end
output.MaxAmp_n = length([output.MaxAmp]);
output.HalfMax_n = length([output.HalfMax]);
output.MaxAmp_SEM = nansem([output.MaxAmp]);
output.HalfMax_SEM = nansem([output.HalfMax]);
end