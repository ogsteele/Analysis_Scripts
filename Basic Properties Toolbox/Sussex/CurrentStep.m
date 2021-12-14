[file,path] = uigetfile('*.ma');
S = ephysIO(fullfile(path,file));
% preallocate pks, locs, w, p
warning('off','signal:findpeaks:largeMinPeakHeight');
Time = S.array(:,1);
Waves = S.array(:,2:end);
for i = 2:size(Waves,2)
    % pks = value of peak
    % locs = index of peak
    % w = width
    % p = prominence
    [pks{i},locs{i},w{i},p{i}] = findpeaks(S.array(:,i),...
        'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',0.02*10^4,...
        'WidthReference','halfheight');
    numSpikes(i) = size(pks{i},1);
end
% plot whole of current step protocol
figure; plot(Time,Waves*1000,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');

% determine rheobase as the first amount of current to induce APs in the
% first 25% of the current step rather than the first current step value to
% elcit any AP at all. 

% First AP of interest values
C = locs;
idx = ~cellfun('isempty',C);
out = zeros(size(C));
out(idx) = cellfun(@(v)v(1),C(idx));
logicalIndexes =  out < 27500 & out > 1;
wavenum_first = (size(locs,2)-(sum(logicalIndexes))); % wave that the first action potential following stimulus injection occurs on
figure; plot(Time,Waves(:,wavenum_first),'color','red'); hold on; plot(Time,Waves(:,11),'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
pA = [-100:10:190]';
lgd = legend(char(string(pA(wavenum_first))),char(string(pA(11))),'linewidth',1);
title(lgd,'Current (pA)')

%% Ih Sag Values
% measured after a depolarising current injection of -100 pA from -65 mV

% Plot Ih Sag Waves
figure; plot(Time,Waves(:,11)*1000,'color','black'); hold on; plot(Time,Waves(:,1)*1000,'color','red')
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
[y,x] = min(Waves(:,1))
x_Time = x*Time(2);
SS_Value = mean(Waves(SS_start:SS_end,1)) + 0.065;
Ih_Sag_Amp = ((SS_Value - (y+0.065)))*1000; % Sag amplitude in mV
Ih_Sag_Percentage = (1-(SS_Value/y))*100; % Sag percentage

hold on; yline((SS_Value-0.065)*1000,'--r'); yline((y*1000),'--r');

%% AP analysis

%% Input Resistance

% plot waves used for input resistance calculation
figure; plot(Time,Waves(:,1)*1000,'color','black'); hold on; plot(Time,Waves(:,2)*1000,'color','red')
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
I = 10e-12; % I (Amps)
R = deltaV / I; % R (Ohms)
IR = R / 1e6; % R (MegaOhms)
txt_1 = ['\bf Input Resistance: '];
txt_2 = [num2str(IR) ' M\Omega \rightarrow'];
hold on; text(0.7,(mean(Waves(IR_start:IR_end,2))*1000) + 5,txt_1)
hold on; text(0.8,(mean(Waves(IR_start:IR_end,2))*1000) + 3,txt_2)