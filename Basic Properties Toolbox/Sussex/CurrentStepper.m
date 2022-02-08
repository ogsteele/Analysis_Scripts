function [output] = CurrentStepper 
%% To Do
% consider finding more robust way of detecting the decay kinetics

%% Code

% user should select the first 'Clamp1.ma' recording in the sequence

% output is a matlab structure containing the following;
% genotype - APOE3 or APOE4
% ketamine exposure - Yes/No
% filepath - filepath of the first wave
% steps - current amplitude of each step in pA
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

% load file of interest
[file,path] = uigetfile('*.ma');
S = ephysIO(fullfile(path,file));
Time = S.array(:,1);
Waves = S.array(:,2:end);

% preallocate pks, locs, w, p, numSpikes and mute error here
warning('off','signal:findpeaks:largeMinPeakHeight');
pks = cell(size(Waves,2),1);
locs = cell(size(Waves,2),1);
w = cell(size(Waves,2),1);
p = cell(size(Waves,2),1);
numSpikes = zeros(size(Waves,2),1);

% findpeaks to determine number and location of AP's per wave
for i = 1:size(Waves,2)
    % pks = value of peak
    % locs = index of peak
    % w = width
    % p = prominence
    [pks{i},locs{i},w{i},p{i}] = findpeaks(Waves(:,i),...
        'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',0.02*10^4,...
        'WidthReference','halfheight');
    numSpikes(i) = size(pks{i},1);
end

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
pA_waveform(0.01*Hz:0.06*Hz,:) = -60; % initial test pulse
for i = 1:size(pA,1)
    pA_waveform(0.5*Hz:1.5*Hz,i) = pA(i); 
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
out = zeros(size(C));
out(idx) = cellfun(@(v)v(1),C(idx));
logicalIndexes =  out< 27500 & out > 1;
wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
subplot(7,4,[2,6,10]);
plot(Time,Waves(:,wavenum_first),'color','red'); hold on; plot(Time,Waves(:,11),'color','black')
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
for i = 1:11
    [y(i),x(i)] = min(Waves(20000:30000,i));
    x_Time = x(i)*Time(2);
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
AP_Window = Waves(out(wavenum_first)-500:out(wavenum_first)+500,wavenum_first)*1000;

% Overshoot in mV
[Overshoot,ind_o] = max(AP_Window);
% Afterhyperpolarisation in mV
[Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% Baseline & Amplitude
Base = mean(AP_Window(1:350));
Amplitude = abs(Base - Overshoot); % in mV

% plot action potential waveform
subplot(7,4,[20,24,28]); plot(AP_Window, gradient(AP_Window),'linewidth',2,'color','black')
box off; title('Action Potential Waveform')
set(gca,'linewidth',2); set(gcf,'color','white'); xlabel('Membrane Potential (mV)'); ylabel('dV/dt');

% Action potential halfwidth
% Halfwidth in ms
subplot(7,4,[4,8,12]);
[AP_pk,AP_l,Halfwidth,AP_pr] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
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
[~,thresh_ind,~,~] = findpeaks(gradient(gradient(gradient(AP_Window))),...
    'MinPeakProminence',0.1, ...
    'MinPeakWidth',3, ...
    'NPeaks',1);
hold on; plot(gradient(gradient(gradient(AP_Window)))*100-20);
Threshold = AP_Window(thresh_ind); % in mV
hold on; plot(thresh_ind,Threshold-Base,'Or') % plot threshold
hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after

% Depolarisation Rate
% between 25 % and 75 % of the rise phase
rise_25 = round(0.25*(ind_o - thresh_ind)) + thresh_ind;
rise_75 = round(0.75*(ind_o - thresh_ind)) + thresh_ind;
Rise = mean(gradient(AP_Window(rise_25:rise_75)));
hold on; plot((rise_25:rise_75),AP_Window(rise_25:rise_75)-Base,'linewidth',2)

% Repolarisation Rate
% between 10 % and 60 % of the falling phase
fall_25 = round(0.1*(ind_o+ind_a - ind_o)) + ind_o;
fall_75 = round(0.6*(ind_o+ind_a - ind_o)) + ind_o;
Fall = mean(gradient(AP_Window(fall_25:fall_75)));
hold on; plot((fall_25:fall_75),AP_Window(fall_25:fall_75)-Base,'linewidth',2,'color','red')

legend('trace','peak', ...
    'amplitude','halfwidth', ...
    'border','dvdt_3', ...
    'threshold','afterhyp', ...
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
txt_1 = ['\bf Input Resistance: '];
txt_2 = [num2str(IR) ' M\Omega \rightarrow'];
hold on; text(0.7,(mean(Waves(IR_start:IR_end,2))*1000) + 5,txt_1)
hold on; text(0.8,(mean(Waves(IR_start:IR_end,2))*1000) + 3,txt_2)

%% Bridge Balance adjustments
% apply necessary bridge balance adjustments if the user requested offline

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


%% Output
% What Genotype was the animal?
dlgTitle    = 'Genotpye';
dlgQuestion = 'What Genotype was this animal?';
genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4','APOE3');

% Was the recording exposed to ketamine?
dlgTitle    = 'Ketamine';
dlgQuestion = 'Was this recorded in the presence of ketamine?';
ketamine = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');

% create output structure
output.genotype = genotype;
output.ketamine = ketamine;
output.filepath = path;
output.steps = pA;
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

% navigate to root dir
cd(path)
cd ..\..

% save output
outname = split(strtrim(path),filesep);
outname = char(string(outname(end-2)));
saveas(fh,[outname,'.fig']);
save([outname,'.mat'],'output')
end