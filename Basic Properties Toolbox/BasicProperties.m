function BasicProperties(capacitance,sweep)
% Author: O.G Steele, 15.10.17
% Description
% Writes to Excel spreadsheet in line with the Cardiff Basic Properties
% Analaysis SOP

% Input Variables
    % capacitance = capacitance in pF as taken from properties of the
    %               ActInact Protocol

    % sweep       = sweep number of the current step protocol containing
    %               the first action potential after injection of current
    
% NOTE 
    % enter sweep number above 21 if no action potential for analysis
%% Variable Code
% Input Rig #
Rig = 5;
    xlswrite('BasicProperties.xlsx',Rig,'Sheet1','A7');
% Input name of Experiment/Project
Experiment = 'RC9 -ACT -X +P';
    xlswrite('BasicProperties.xlsx',cellstr(Experiment),'Sheet1','A4');
% Input DIV
DIV = 27;
    xlswrite('BasicProperties.xlsx',DIV,'Sheet1','A5');
% Input Experimentor
Experimentor = 'O.G. Steele';
    xlswrite('BasicProperties.xlsx',cellstr(Experimentor),'Sheet1','A1');

%% Code
abf = dir('*.abf'); 
numfiles = length(abf);
mydata = cell(1, numfiles);

for k = 1:numfiles 
  mydata{k} = abfload(abf(k).name); 
end
GapFree = mydata{1,1};
CurrentStep = mydata{1,3};
ActInact = mydata{1,4};

abf1 = abf(1);
    abf1name = getfield(abf1,'name');
        xlswrite('BasicProperties.xlsx',cellstr(abf1name),'Sheet1','A6');
    abf1date = getfield(abf1,'date');
        xlswrite('BasicProperties.xlsx',cellstr(abf1date),'Sheet1','A2');
    abf1folder = getfield(abf1,'folder');
        xlswrite('BasicProperties.xlsx',cellstr(abf1folder),'Sheet1','A3');
abf2 = abf(2);
abf3 = abf(3);
    abf3name = getfield(abf3,'name');
        xlswrite('BasicProperties.xlsx',cellstr(abf3name),'Sheet1','A28');
abf4 = abf(4);

% Write in filename1, capacitance
xlswrite('BasicProperties.xlsx',capacitance,'Sheet1','A10');

% Resting Membrane Potential (Vm)(measured in mV) Data
% Where Vm is mean of the first 10ms of gapfree recording
Vm = mean(GapFree(1:200000));
    xlswrite('BasicProperties.xlsx',Vm,'Sheet1','A8');
    
% Spontaneous Action Potential Code
%   - Spontaneous defined as above 0mV
%   - Attempted defined as above Vm + 10mV)
%   - Quiet defined as everything else
if GapFree > 0
    xlswrite('BasicProperties.xlsx','3','Sheet1','A17');
    xlswrite('BasicProperties.xlsx','3','Sheet1','A18');
elseif GapFree > (Vm + 5)
     xlswrite('BasicProperties.xlsx','2','Sheet1','A16');
     xlswrite('BasicProperties.xlsx','2','Sheet1','A18');
else
    xlswrite('BasicProperties.xlsx','1','Sheet1','A14');
    xlswrite('BasicProperties.xlsx','1','Sheet1','A18');
end
    
% Input Resistance (IR)(measured in MOhm) Data
% Difference of the mean of the last 10ms prior to restoration to -70mV or
% of quiet section where traces are stable
time = (1:103224);
plot(time,CurrentStep(:,1:2))
hold on
    title('IR Calculator Point 1 = start, Point 2 = end')
    hold off
[x] = getpts;
    a = round(x(2));
    b = round(x(1));
IR = (mean(CurrentStep(b:a,:,1) - CurrentStep(b:a,:,2)))*-100;
    xlswrite('BasicProperties.xlsx',IR,'Sheet1','A9');  

% Sodium Activation Data
% CONSTANTS
pA = linspace(-120,60,37)';
DF = pA-66.68;
% Prerequisites
DavidActInactPro = reshape(ActInact,10322,37);
% Calculations
NaAVal = DavidActInactPro(1167:1505,:,:);
NaAPeaks = (min(NaAVal))';
    xlswrite('BasicProperties.xlsx',NaAPeaks,'Sheet1','A86');
NaADensity = NaAPeaks/capacitance;
    xlswrite('BasicProperties.xlsx',NaADensity,'Sheet1','A124');
MinNaADensity = min(NaADensity);
     xlswrite('BasicProperties.xlsx',MinNaADensity,'Sheet1','A82');
NaAConductance = NaADensity./DF;
    xlswrite('BasicProperties.xlsx',NaAConductance,'Sheet1','A163');
NaANormConductance = NaAConductance./max(NaAConductance);
    xlswrite('BasicProperties.xlsx',NaANormConductance,'Sheet1','A203')

% Sodium Inactivation Data
NaIVal = DavidActInactPro(5000:6000,:,:);
NaIPeaks = (min(NaIVal))';
    xlswrite('BasicProperties.xlsx',NaIPeaks,'Sheet1','A243');
NaIConductance = NaIPeaks./(-66.68);
    xlswrite('BasicProperties.xlsx',NaIConductance,'Sheet1','A282');
NaINormConductance = NaIConductance/max(NaIConductance);
    xlswrite('BasicProperties.xlsx',NaINormConductance,'Sheet1','A322');

% Potassium Current Data
KVal = DavidActInactPro(4160:5160,:,:);
KMean = (mean(KVal))';
    xlswrite('BasicProperties.xlsx',KMean,'Sheet1','A361');
KDensity = KMean/capacitance;
    xlswrite('BasicProperties.xlsx',KDensity,'Sheet1','A400');
MaxKDensity = max(KDensity);
    xlswrite('BasicProperties.xlsx',MaxKDensity,'Sheet1','A83');

% Spike Analysis
if sweep > 21
    NaN = 0;
    xlswrite('BasicProperties.xlsx',NaN,'Sheet1','A28:A36');
else
    SpikeSweepWindow = CurrentStep(1:10000,:,sweep);
    time = (1:10000)';
    plot(time,SpikeSweepWindow);
        hold on
        title('Define the AP, where P1 = pre, P2 = peak, P3 = post')
        hold off
    [x] = getpts;
       pre = round(x(1));
       peak = round(x(2));
       post = round(x(3));
   
    PeakmV = max(SpikeSweepWindow(pre:post,:,:));
       xlswrite('BasicProperties.xlsx',PeakmV,'Sheet1','A31');
    AftHyp = min(SpikeSweepWindow(peak:post,:,:));
       xlswrite('BasicProperties.xlsx',AftHyp,'Sheet1','A32');
    SpikeAmp = PeakmV - AftHyp;
       xlswrite('BasicProperties.xlsx',SpikeAmp,'Sheet1','A33'); 
    Rise = (max(gradient(smooth(SpikeSweepWindow(pre:peak)))))*200;
       xlswrite('BasicProperties.xlsx',Rise,'Sheet1','A34');
    Fall = (min(gradient(smooth(SpikeSweepWindow(peak:post)))))*200;
       xlswrite('BasicProperties.xlsx',Fall,'Sheet1','A35');
    HalfWHeight = 0 - (SpikeAmp/2);
    plot (SpikeSweepWindow(pre:post));
       hold on
        title('Half Width where P1 = rising intercept, P2 = falling intercept')
        hline = refline(0,HalfWHeight);
       hold off
    [x] = getpts;
        HalfStart = round(x(1));
        HalfEnd = round(x(2));
    HalfWidth = (HalfEnd - HalfStart) / 50;
        xlswrite('BasicProperties.xlsx',HalfWidth,'Sheet1','A36');
    Smo = smooth(SpikeSweepWindow);
    time = (1:9000)';
    CS10win = Smo(1:peak);
    timewin = (1:peak);
    TDiff = smooth(diff(diff(diff(CS10win))));
    TimeDiff = (1:(peak-3))';
    plot(timewin,CS10win,TimeDiff,TDiff);
        hold on
         title('Threshold Potential where P1 = pre threshold, P2 = post threshold')
        hold off
    [x,y] = getpts;
        pre = round(x(1));
        post = round(x(2));
        close Figure 1
    TDiff2 = TDiff(pre:post);
    [a,b] = max(TDiff2);
    CS10win2 = Smo(pre:post);
    Thresh = CS10win2(b);
        xlswrite('BasicProperties.xlsx',Thresh,'Sheet1','A30');
        xlswrite('BasicProperties.xlsx',sweep,'Sheet1','A29');  
end
end
