function [output] = IStep_K
%% About
%@OGSteele, 2022

% user should select the first 'Clamp1.ma' recording in the sequence

% output is a matlab structure containing the following;
% genotype - APOE3 or APOE4
% ketamine exposure - Yes/No
% filepath - filepath of the first wave
% ID - slice_ID on plate
% steps - current amplitude of each step in pA
% waveform - current step protocol waveform
% time - time in s
% waves - waves in V
% ephysIO - ephysIO output format, with associated metadata if needed
% numSpikes - number of spikes per wave
% Rh - Rheobase in pA
% sag_mV - Ih Sag Amplitude in mV
% sag_ratio - Ih Sag : steady state rati
% peak - Overshoot in mV
% afterhyp - Afterhyperpolarisation value in mV
% amp - Action potential amplitude in mV
% thresh - Threshold potential in mV
% half - Halfwidth in ms
% rise - Depolarisation rate in mV/s
% fall - Repolarisation rate in mV/s
% IR - Input Resistance in MOhm
% Rs_Init - Access resistance approximated from initial step  in MOhm
% Offline_BB - Vm adjustments for offline bridge balance in V
% Online_BB_performed = Yes/No 
% Offline_BB_performed = Yes/No 

% bridge balance corrected values
    % IR
    % thresh
    % afterhyp
    % peak

%% Code

% load file of interest
[file,path] = uigetfile('*.ma');
S = ephysIO(fullfile(path,file));
Time = S.array(:,1);
Waves = S.array(:,2:end);

% cd to path
splitPath = split(path,filesep);
newPath = char(string(join(splitPath(1:end-3),"\")));
cd(newPath)

% check to see if I_step has already been run
try
    list = dir('**\*.mat'); % ending in .mat
    pat = 'CC_IV_1_full_';
    for i = 1:size(list,1)
        list(i).isfile = startsWith(list(i).name,pat); % is it one of ours?
    end
    fileFlags = [list.isfile];
    list = list(fileFlags); % save only the ones we're interested in
    istep_run = 1;
    ('IStep run detected, running with extracted settings as defaults')
    load(list.name)
catch
    warning('Assuming IStep has not been run, running IStep normally');
    istep_run = 0;
end

% run analysis or not?
if istep_run == 0
    % plot figure
    fh = figure();
    plot(Time,Waves*1000,'color','black')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Membrane Potential (mV)'); xlabel('Time (s)')
    title('Current Step Waveform')
    
    % ask the user 
    dlgTitle    = 'Run Analysis';
    dlgQuestion = 'Would you like to run the Current Step Analysis?';
    run = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
    
    close(fh) % close figure
elseif istep_run == 1
    run = "Yes";
end
    % establish if statement for protocol
    if run == "Yes"
        disp('Performing Analysis, please wait ...')

% preallocate pks, locs, w, p, numSpikes and mute error here
warning('off','signal:findpeaks:largeMinPeakHeight');
pks = cell(size(Waves,2),1);
locs = cell(size(Waves,2),1);
difflocs = cell(size(Waves,2),1); % difference between loc in dp
normlocs = cell(size(Waves,2),1); % difference between loc in dp
w = cell(size(Waves,2),1);
p = cell(size(Waves,2),1);
numSpikes = zeros(size(Waves,2),1);
nlocs = cell(size(Waves,2),1); % interval index


% findpeaks to determine number and location of AP's per wave
for i = 1:size(Waves,2)
    % pks = value of peak
    % locs = index of peak
    % w = width
    % p = prominence
    [pks{i},locs{i},w{i},p{i}] = findpeaks(Waves(:,i),...
        'MinPeakHeight',0,...
        'MinPeakProminence',0.01,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',0.02*10^4,...
        'WidthReference','halfheight');
    numSpikes(i) = size(pks{i},1);
    difflocs{i} = diff(locs{i}); % difference between loc in dp
    normlocs{i} = normalize(difflocs{i},'scale','first'); % difflocs norm to first value
    nlocs{i} = 1:size(difflocs{i}); % interval index
end


% interval_index = [];
% normalised_interval = [];
% for i = 1:size(nlocs,1)
%     interval_index = [interval_index,nlocs{i}];
%     normalised_interval = [normalised_interval,normlocs{i}'];
% end
% [p,S] = polyfit(interval_index,normalised_interval,1); % coefficients of first degree polynomial
% [y_fit,delta] = polyval(p,interval_index,S); % fit the polynomial and error
% figure; plot(interval_index,normalised_interval,'bo')
% hold on
% plot(interval_index,y_fit,'r-')
% plot(interval_index,y_fit+2*delta,'m--',interval_index,y_fit-2*delta,'m--')
% title('Linear Fit of Data with 95% Prediction Interval')
% legend('Data','Linear Fit','95% Prediction Interval')
% box off; set(gcf,'color','white'); set(gca,'linewidth',2)
% xlabel('AP Interval #'); ylabel('Normalised Inter Event Interval')


% plot whole of current step protocol
fh = figure();
fh.WindowState = 'maximized'; subplot(7,4,[1,5]); plot(Time,Waves*1000,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
ylabel('Membrane Potential (mV)');
title('Current Step Waveform')
ax = gca; xax = ax.XAxis; set(xax,'visible','off')


% plot the waveform of interest
wavelength = Time(end);
Hz = size(S.array,1)/wavelength;
pA = linspace(-200,400,size(Waves,2))';
pA_waveform = zeros(size(Waves,1),size(Waves,2)); 
pA_waveform(round(0.01*Hz):round(0.06*Hz),:) = -60; % initial test pulse
for i = 1:size(pA,1)
    pA_waveform(round(0.5*Hz):round(1.5*Hz),i) = pA(i); 
end
subplot(7,4,9); plot(Time,pA_waveform(:,:),'linewidth',1,'color','black')
box off; set(gca,'linewidth',2); set(gcf,'color','white');
xlabel('Time (s)'); ylabel('Command(pA)'); ylim([-250 300])
ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')

% determine rheobase as the first amount of current to induce APs in the
% first 25% of the current step rather than the first current step value to
% elcit any AP at all. 

% First AP of interest values
C = locs;
idx = ~cellfun('isempty',C);
outs = zeros(size(C));
outs(idx) = cellfun(@(v)v(1),C(idx));
logicalIndexes =  outs< 27500 & outs > 1;
wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
subplot(7,4,[2,6,10]);
plot(Time,Waves(:,wavenum_first)*1000,'color','red'); hold on; plot(Time,Waves(:,11)*1000,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
% pA = [-100:10:190]'; % legacy numbers
pA = linspace(-200,400,size(Waves,2))';

lgd = legend(char(string(pA(wavenum_first))),char(string(pA(11))),'linewidth',1);
title(lgd,'Current (pA)')
title('Exemplary Waves')

%% Membrane Excitability and Rheobase
% determine rheobase as the first amount of current to induce APs in the
% first 25% of the current step rather than the first current step value to
% elcit any AP at all. 

Rh = pA(wavenum_first); % Rheobase in pA

% plot membrane excitability and Rheobase
subplot(7,4,[3,7,11])
plot(pA,numSpikes,'color','red','linewidth',3); box off; set(gcf,'color','white'); set(gca,'linewidth',2);
title('Membrane Excitability'); xlabel('Current Step (pA)'); ylabel('Number of Action Potentials')
hold on; xline(Rh,'--','linewidth',1.5); 
txt_Rh_1 = ['\bf Rheobase:  ' num2str(Rh) ' pA \rightarrow'];
text(Rh-100,8,txt_Rh_1);

%% Ih Sag Values
% measured after a depolarising current injection of -100 pA from -65 mV
% Calculates both Amplitude (relative to steady state) and ratio (to steady state)

% Plot Ih Sag Waves
subplot(7,4,[17,21,25]); plot(Time,Waves(:,11)*1000,'color','black'); hold on; plot(Time,Waves(:,1)*1000,'color','red')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
lgd = legend(char(string(pA(11))),char(string(pA(1))),...
    'linewidth',1,...
    'location','southeast',...
    'AutoUpdate','off');
title(lgd,'Current (pA)')
title('Ih Sag Calculation')
hold on; plot(Time, Waves(:,2:10)*1000, 'color',[0.8,0.8,0.8])
hold on; plot(Time,Waves(:,11)*1000,'color','black'); hold on; plot(Time,Waves(:,1)*1000,'color','red')
%hold on; xline(1.5,'--'); xline(1.3,'--')

% Calculate Ih Sag
SS_start = 1.3/Time(2); 
SS_end = 1.5/Time(2);
% preallocate for loop variables
y = zeros(1,11);
x = zeros(1,11);
SS_Value = zeros(1,11);
Ih_Sag_Amp = zeros(1,11);
Ih_Sag_Percentage = zeros(1,11);
for i = 1:11
    [y(i),x(i)] = min(Waves(20000:30000,i));
    SS_Value(i) = mean(Waves(SS_start:SS_end,i)) + 0.065;
    Ih_Sag_Amp(i) = ((SS_Value(i) - (y(i)+0.065)))*1000; % Sag amplitude in mV
    Ih_Sag_Percentage(i) = ((SS_Value(i)/y(i)))*100; % Sag percentage
end 
hold on; yline((SS_Value(1)-0.065)*1000,'--r'); yline((y(1)*1000),'--r');

subplot(7,4,[18,22,26]);
plot(pA(1:11),Ih_Sag_Percentage,'color','black','linewidth',3); box off; title('Ih Sag')
set(gca,'linewidth',2); set(gcf,'color','white'); xlabel('Current Step (pA)'); ylabel('Ih Sag - Steady State Ratio (%)');
%% AP analysis

% plot action potential
AP_Window = Waves(outs(wavenum_first)-500:outs(wavenum_first)+500,wavenum_first)*1000;

% Overshoot in mV
[Overshoot,ind_o] = max(AP_Window);
% Afterhyperpolarisation in mV
[Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% Baseline
Base = mean(AP_Window(1:350));

% plot action potential waveform
subplot(7,4,[20,24,28]); plot(AP_Window, gradient(AP_Window),'linewidth',2,'color','black')
box off; title('Action Potential Waveform')
set(gca,'linewidth',2); set(gcf,'color','white'); xlabel('Membrane Potential (mV)'); ylabel('dV/dt');

% Action potential halfwidth
% Halfwidth in ms
subplot(7,4,[4,8,12]);
[~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
set(gca,'linewidth',2); set(gcf,'color','white'); title('Action Potential Analysis');
ylim([-40 120])
Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms
%deal with the legend later

% Action Potential Threshold
% as defined by the first peak of the third derivative
% (or second derivative of dv/dt as per Serkeli et al., 2004)
% what I actually prefer is when the dv/dt surpasses 4mV/ms (dv/dt > 1)
N = diff(AP_Window(1:ind_o-20)); % period start to just before the peak
[closestValue, closestIndex] = min(abs(N - 1.'));
hold on; plot(closestIndex,AP_Window(closestIndex)-Base,'or')
ind_t = closestIndex; % threshold index
%hold on; plot(gradient(gradient(gradient(AP_Window)))*100-20);
% plot the differntial below the trace
hold on; plot(diff(AP_Window)*1-20);
hold on; yline(closestValue-20) % threshold (dv/dt)
Threshold = AP_Window(ind_t); % in mV

% PennThresh
% peak of the first derivative

%hold on; plot(thresh_ind,Threshold-Base,'Or') % plot threshold
hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after

% recalculate (and overwrite) the amplitude now you have a threshold
Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV

% Depolarisation Rate
% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
rise_20_ind = closestIndex_20 + ind_t;
rise_80_ind = closestIndex_80 + ind_t;
Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)

% Repolarisation Rate

% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
fall_20_ind = closestIndex_80 + ind_o;
fall_80_ind = closestIndex_20 + ind_o;
Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)

legend('trace','peak', ...
    'amplitude','halfwidth', ...
    'border','threshold', ...
    'dv/dt','dv/dt > 1',...
    'hyperpolar.',...
    'rise','fall',...
    'Location','northeast')
%% Input Resistance
% Calculates input resistance in MOhm from the difference in voltage, divided by
% the current applied, to the the steady state potentials on the last two
% waves

% plot waves used for input resistance calculation
subplot(7,4,[19,23,27]); plot(Time,Waves(:,1)*1000,'color','black'); hold on; plot(Time,Waves(:,2)*1000,'color','red')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
lgd = legend(char(string(pA(1))),char(string(pA(2))),...
    'linewidth',1,...
    'location','southeast',...
    'AutoUpdate','off');
title(lgd,'Current (pA)')
title('Input Resistance Calculation')
%show lines for region of IR determination
hold on; xline(1.5,'--'); xline(1.3,'--')

% Calculate Input Resistance in MegaOhms
IR_start = 1.3/Time(2); 
IR_end = 1.5/Time(2);
deltaV = abs(mean(Waves(IR_start:IR_end,1)) - mean(Waves(IR_start:IR_end,2))); % Delta_Voltage (Volts)
I = 20e-12; % I (Amps)
R = deltaV / I; % R (Ohms)
IR = R / 1e6; % R (MegaOhms)
txt_1 = '\bf Input Resistance: ';
txt_2 = [num2str(IR) ' M\Omega \rightarrow'];
hold on; text(0.7,(mean(Waves(IR_start:IR_end,2))*1000) + 5,txt_1)
hold on; text(0.8,(mean(Waves(IR_start:IR_end,2))*1000) + 3,txt_2)

%% Bridge Balance adjustments
% apply necessary bridge balance adjustments if the user requested offline

if istep_run == 0
% plot the initial current step so user can see if recording was balanced
t = figure; plot(Time(1:151),Waves(300:450,1)); box off; xlabel('Time (s)');
title('Initial Current Step [Zoomed]'); set(gca,'linewidth',2); 
set(gcf,'color','white'); ylabel('Membrane Potential (mV)'); 

% was this recording bridge balanced appropriately?
dlgTitle    = 'Bridge Balance';
dlgQuestion = 'Was this recording appropriately bridge balanced?';
balanced = questdlg(dlgQuestion,dlgTitle,'Yes','No','Yes');

close(t) % closes the test fig as it's not interesting anymore

% logical fork following balancing
if balanced == "Yes"
    disp('Recording appropriately balanced, no further action required')
    Vm_adjust = NaN(size(Waves,2),1);
    Rs_Init = NaN;
    offline_BB_performed = "No";
elseif balanced == "No"
    disp('Recording not appropriately balanced, performing offline bridge balance')
    offline_BB_performed = "Yes";
    % calculate Rs values from the initial current step
    Istep = -60; % pA
    % find the triple diff peak
    [~, ind] = findpeaks((gradient(gradient(Waves(:,1))))*-1,...
        'minPeakProminence',2.5e-4,'NPeaks',1);
    % preallocate
    Vm_1 = zeros(size(Waves,2),1);
    Vm_2 = zeros(size(Waves,2),1);
    delta_Vm = zeros(size(Waves,2),1);
    Vm_adjust = zeros(size(Waves,2),1);
    for i = 1:size(Waves,2)
        Vm_1(i) = mean(Waves(1:ind,i)); % in V
        Vm_2(i) = mean(Waves(405:410,i)); % in V
        delta_Vm(i) = abs(Vm_1(i)-Vm_2(i))*1000; % in mV
    end
    Rs_Init = abs(median(delta_Vm/Istep)*1000); % in MOhm
    for i = 1:size(Waves,2)
        Vm_adjust(i) = ((Rs_Init*1e6)*(pA(i)*1e-12)); % in V (initial)
    end
    
    % plot the required adjustment
    figure; plot(pA, Vm_adjust*1000,'-o'); box off; set(gca,'linewidth',2); 
    set(gcf,'color','white'); xlabel('pA'); ylabel('Adjustment req. (mV)')
    title('Offline Bridge Balance Adjustment Required')
    
    % balance the data here
    % IR adjustment
    deltaV = abs((mean(Waves(IR_start:IR_end,1))-Vm_adjust(1)) ...
        - (mean(Waves(IR_start:IR_end,2))-Vm_adjust(2))); % Delta_Voltage (Volts)
    I = abs(pA(1)*1e-12 - pA(2)*1e-12); % I (Amps)
    R = deltaV / I; % R (Ohms)
    IR = R / 1e6; % R (MegaOhms)
    % thresh
    Threshold = Threshold-(Vm_adjust(wavenum_first)*1000);
    % afterhyp
    Afterhyperpolarisation = Afterhyperpolarisation-(Vm_adjust(wavenum_first)*1000);
    % peak
    Overshoot = Overshoot-(Vm_adjust(wavenum_first)*1000);

end
    elseif istep_run == 1
    balanced = output.Online_BB_performed;
    Rs_Init = output.Rs_Init;
    Vm_adjust = output.Offline_BB;
    IR = output.IR;
    Threshold = output.thresh;
    Afterhyperpolarisation = output.afterhyp;
    Overshoot = output.peak;
    offline_BB_performed = output.Offline_BB_performed;
end
%% Output
% What Genotype was the animal?
dlgTitle    = 'Genotpye';
dlgQuestion = 'What Genotype was this animal?';
if istep_run == 0
    genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4','APOE3');
elseif istep_run == 1
    genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4',(output.genotype));  
end
% Was the recording exposed to ketamine?
dlgTitle    = 'Ketamine';
dlgQuestion = 'Was this recorded in the presence of ketamine?';
if istep_run == 0
    ketamine = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
elseif istep_run == 1
    ketamine = questdlg(dlgQuestion,dlgTitle,'Yes','No',(output.ketamine));
end

% What was the slice ID, and  add any notes?
prompt = {'What was the slice ID for this recording?','Would you like to add any notes?'};
dlgtitle = 'Slice ID & notes';
if istep_run == 0
    definput = {'eg. E420211203#1','eg. dodgy input'};
elseif istep_run == 1
    definput = {char(output.ID),char(output.Notes)};
end
dims = [1 40];
slice_id_notes = inputdlg(prompt,dlgtitle,dims,definput);


% Would you like to add any notes? 
% prompt = {'Would you like to add any notes'};
% dlgtitle = 'Notes';
% definput = {'eg. dodgy input'};
% dims = [1 40];
% notes = inputdlg(prompt,dlgtitle,dims,definput);

% create output structure
output.genotype = genotype;
output.ketamine = ketamine;
output.filepath = path;
output.ID = slice_id_notes(1);
output.steps = pA;
output.waveform = pA_waveform;
output.time = Time;
output.waves = Waves;
output.ephysIO = S;
output.numSpikes = numSpikes;
output.Rh = Rh;
output.sag_mV = Ih_Sag_Amp;
output.sag_ratio = Ih_Sag_Percentage;
output.peak = Overshoot;
output.afterhyp = Afterhyperpolarisation;
output.amp = Amplitude;
output.thresh = Threshold;
output.half = Halfwidth;
output.rise = Rise;
output.fall = Fall;
output.IR = IR;
output.Rs_Init = Rs_Init;
output.Offline_BB = Vm_adjust;
output.Online_BB_performed = balanced;
output.Offline_BB_performed = offline_BB_performed;
output.Notes = slice_id_notes(2);

% navigate to root dir
cd(path)
cd ..\..

% save output
outname = split(strtrim(path),filesep);
outname = char(string(outname(end-2)));
saveas(fh,[outname,'.fig']);
save([outname,'.mat'],'output')

% return to if loop from the top 
    else
        close(fh) % close figure
        disp('Analysis Aborted, please find a new recording')
    end
end