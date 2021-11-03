%% Rs_Comp_OGS
% Script to load in ephys data, and compensate the series resistance
% changes observed during whole cell patch clamp recordings

% Note, change parameters as relevent

% Dependencies - ensure in path
% ephysIO (https://github.com/acp29/eventer/blob/master/base/ephysIO.m)


% Output:
% .parameters = experimental parameters
% .raw_splits = raw recording split into 10s waves 
% .raw_wcp = raw whole cell properties
% .raw_ephysIO = ephysIO output from raw trace
% .comp_splits = rs comp recording split into 10s waves 
% .comp_wcp = rs comp whole cell properties
% .comp_ephysIO = ephysIO output from rs comp trace

% figures:
% raw wcp values and filtered trace
% compensated wcp values and filtered trace
%% TO DO
% Consider final output structures
% Change units on axis 7 of master subplot
% add units and titles to subplots

% start with a blank workspace 
clear
close all

%% Options
% for speed, set plots to off
plots = false; % true / false (logical)
comp_penn = true; % leave as true (logical)

%% Parameters
% Amplifier / Sampling settings
Param.amp_scalef = 0.5; % default scale factor for amplifier used in V / nA
Param.sample_rate = 20000; % in Hz
Param.split_length = 10; % in seconds, frequency of test pulses
% Test pulse settings
Param.pulse_amp = -0.002; % in V
Param.Vh = -45; % holding voltage in mV
Param.voltage_step = 2; % in mV
Param.pulse_duration = 10; % in ms
Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
% Compensation settings
Param.Vrev = 0; % reversal potential in mV (close to zero for AMPAR/NMDAR)
Param.des_Rs = 8; % desired series resistance for recordings to be compensated to

%% Load in w/ ePhysIO
% Select raw trace to visualise
title_str = "1. Select raw file of recording";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[file,path,~] = uigetfile('*.*','1. Select raw file of recording');
% Display file selection selection
if isequal(file,0)
   disp('User selected Cancel')
   % If user selects cancel here, script will end here.
   return
else
    disp(['User selected ', fullfile(path, file)])
    filename = file;
    % Navigate to directory and load file with ephysIO
    cd(path)
    S = ephysIO(file);
end
% tidy workspace
clear('path','ans')

%% What gain level was the recording performed at
prompt = {'Gain value of the raw recording'};
dlgtitle = 'Input Gain setting';
dims = [1 50];
def = {'100'}; % default gain value in my case is 100
gain = inputdlg(prompt,dlgtitle,dims,def);
gain = str2double(gain); % convert to number

% overwrite the gain value set in parameters
Param.amplifier_gain = gain;

%% What Vm was the recording performed at
prompt = {'Holding potential in mV'};
dlgtitle = 'Holding Potential (mV)';
dims = [1 50];
def = {'-70'}; 
Vh = inputdlg(prompt,dlgtitle,dims,def);
Vh = str2double(Vh); % convert to number

% overwrite the gain value set in parameters
Param.Vh = Vh;
%% Split the data into ten second waves

% convert the data in pA
S.array(:,2) = (S.array(:,2)/(Param.amp_scalef * Param.amplifier_gain)*1000); % in pA
% calculate the length of the recording
length_raw = length(S.array(:,2)); % whole recording in data points
length_seconds = length_raw/Param.sample_rate; % whole recording in seconds
n_splits = floor(length_seconds/Param.split_length); % number of splits, rounded down 
% calculate the number of data points per split
length_raw_split = Param.split_length * Param.sample_rate;

% create list of start and end points determined by the sample rate, length
% of split, length of recording
% round down n_splits to not deal with incomplete sections and preallocate
start(:,1) = zeros(n_splits,1);
finish(:,1) = zeros(n_splits,1);
% define the starting numbers
start(1,1) = 1;
finish(1,1) = 1 + length_raw_split-1;
for i = 2:n_splits
    start(i,1) = finish(i-1,1) + 1;
    finish(i,1) = start(i,1) + length_raw_split-1;
end

% create a matrix of the split recording and preallocate
splits = zeros(length_raw_split,n_splits);
for s = 1:n_splits
    splits(:,s) = S.array(start(s):finish(s),2);
end

% tidy up workspace
clear('start','finish','i','s','length_seconds','length_raw','length_raw_split','n_splits')

% get the start time of the test pulse
figure; plot(splits(:,1))
title("2. Zoom into the test pulse BEFORE clicking enter to select a single start point")
pause
title_str = "2. Zoom into the test pulse BEFORE clicking enter to select a single start point";
if ~ispc; menu(title_str,'OK'); end
clear('title_str')
[x,y] = ginput(1);
close

% create pulse paramaters from ginput selection
Param.pulse_start = round(x);
Param.pulse_end = Param.pulse_start + Param.pulse_points;                                     % note, this will only cover the downward transiet, not the upward transient
Param.pulse_window = Param.pulse_start:Param.pulse_end;

% plot first and last window used to calculate wcp parameters
%figure; plot(splits(Param.pulse_window,1))
%hold on
%plot(splits(Param.pulse_window,size(splits,2)))
%legend("first pulse","last pulse")
%title("first and last raw test pulses overlaid")

%% Median filter the recording minus the test pulse
poles = 9;
perf_noise = menu('Was there perfusion noise present?','Yes','No');
if perf_noise == 1 
    [splits, nf_splits] = TP_null_MedianF_function(Param,splits,x,poles);
else
end

%% Generate necessary whole cell paramaters
wcp_raw = WCP(splits,Param);

%% Compensation Penn
% adapted from rscomp_Penn.m

if comp_penn == true
    fraction = 1 - (Param.des_Rs ./ wcp_raw.Rs);
    disp("Comp_penn enabled (to disable set comp_penn = false)")
    % preallocate for speed
    corr_p_splits = zeros(size(splits,1),size(splits,2));
    for i = 1:size(splits,2)
      % y must a vector of the current trace (in pA)
      % sampInt is sampling interval (in microseconds)
      % fraction is the amount of compensation (1.0 = 100%)
      % R_s is series resistance (in Mohm)
      % C_m is cell time constant C_m (in ms)
      % V_hold is holding potential (in mV)
      % V_rev is reversal potential of the conductance (in mV)
      % eg y = rscomp (y, sampInt, fraction, R_s, C_m, V_hold, V_rev)
      sampInt = 1e6 / Param.sample_rate;
      R_s = wcp_raw.Rs(i);
      C_m = wcp_raw.Cm(i);
      V_hold = Param.Vh;
      corr_p_splits(:,i) = rscomp (splits(:,i), sampInt, fraction(i), R_s, C_m, V_hold);
    end
else
    disp("Comp_penn disabled (to enable set comp_penn = true)")
end

% calculate whole cell properties from the compensated traces
wcp_corr_penn = WCP(corr_p_splits,Param);
% save the compensated trace
s(:,2) = (vertcat(corr_p_splits(:))) * 1e-12;
t = 0:(1/Param.sample_rate):(1/Param.sample_rate)*(length(s));
s(:,1) = t(1:end-1);
a = split(filename,'.');
filename_comp = append(char(a(1)),'_rscomp.phy');
ephysIO(filename_comp,s,S.xunit,S.yunit,S.names,S.notes,'int16')


% output
comp_output.parameters = Param;
comp_output.raw_splits = splits; 
comp_output.raw_wcp = wcp_raw;

S.array(:,2) = S.array(:,2)*1e-12; % convert to amps for saving
comp_output.raw_ephysIO = S;
comp_output.comp_splits = corr_p_splits; 
comp_output.comp_wcp = wcp_corr_penn;

S.array = s; % swap raw array with compensated array
comp_output.comp_ephysIO = S;

% save as name of file in dir
a = split(file,'.');
file = append(char(a(1)),'_wcp.mat');
% save(file,'comp_output')
% if larger than 2GB error will arise, use this line instead
save(file, 'comp_output', '-v7.3')
% NOTE TO CHANGE BACK AFTER

% tidy up
clear sampInt R_s C_m V_hold s t comp_penn filename filename_comp a 
%% Plots
if plots == true
    disp("Plots enabled (to disable set plots = false)")
    [YF, median_YF, x, x2] = TestPulse_rm(splits,Param);
    % holding current (Ih) plot
    figure
    set(gcf,'color','w');
    ax1 = subplot(3,3,1);
    plot(ax1, wcp_raw.Ih)
    title('Ih in pA')
    hold on; plot(ax1, movmean(wcp_raw.Ih,7)); hold off
    % series resistance (Rs) plot
    ax2 = subplot(3,3,2);
    plot(wcp_raw.Rs)
    title('Rs in MOhm')
    hold on; plot(movmean(wcp_raw.Rs,7)); hold off
    % input resistance (Rin) plot
    ax3 = subplot(3,3,3);
    plot(wcp_raw.Rin)
    title('Rin in MOhm')
    hold on; plot(movmean(wcp_raw.Rin,7)); hold off
    % memrane resistance (Rm) plot
    ax4 = subplot(3,3,4);
    plot(wcp_raw.Rm)
    title('Rm in MOhm')
    hold on; plot(movmean(wcp_raw.Rm,7)); hold off
    % capacitace (Cm) plot
    ax5 = subplot(3,3,5);
    plot(wcp_raw.Cm)
    title('Cm in pF')
    hold on; plot(movmean(wcp_raw.Cm,7)); hold off
    % steady state current (Iss) plot
    ax6 = subplot(3,3,6);
    plot(wcp_raw.Iss)
    title('Iss in pA')
    hold on; plot(movmean(wcp_raw.Iss,7)); hold off
    % plot without test pulse
    ax7 = subplot(3,3,[7,9]);
    plot(x,YF,'Color',[0 0.4470 0.7410 0.2]) % to make lighter
    title('Recording w/ 1 kHz LPF')
    xlabel('Time (s)')
    ylabel('Current (pA)')
    hold on
    plot(x2,smooth(median_YF),'linewidth',2,'color',[0 0.4470 0.7410])
    hold off % not made lighter
    sgtitle('Raw WCP Values by epoch')
    
    % set axis properties
    ax1.LineWidth = 2;
    ax1.Box = 'off';
    ax2.LineWidth = 2;
    ax2.Box = 'off';
    ax3.LineWidth = 2;
    ax3.Box = 'off';
    ax4.LineWidth = 2;
    ax4.Box = 'off';
    ax5.LineWidth = 2;
    ax5.Box = 'off';
    ax6.LineWidth = 2;
    ax6.Box = 'off';
    ax7.LineWidth = 2;
    ax7.Box = 'off';
    
    % save figure1
    a = split(file,'.');
    file_raw = append(char(a(1)),'_raw.pdf');
    saveas(gcf,file_raw)
    
    % create figure 2
    [YF, median_YF, x, x2] = TestPulse_rm(corr_p_splits,Param);
    % holding current (Ih) plot
    figure
    set(gcf,'color','w');
    ax1 = subplot(3,3,1);
    plot(ax1, wcp_corr_penn.Ih)
    title('Ih in pA')
    hold on; plot(ax1, movmean(wcp_corr_penn.Ih,7)); hold off
    % series resistance (Rs) plot
    ax2 = subplot(3,3,2);
    plot(wcp_corr_penn.Rs)
    title('Rs in MOhm')
    hold on; plot(movmean(wcp_corr_penn.Rs,7)); hold off
    % input resistance (Rin) plot
    ax3 = subplot(3,3,3);
    plot(wcp_corr_penn.Rin)
    title('Rin in MOhm')
    hold on; plot(movmean(wcp_corr_penn.Rin,7)); hold off
    % memrane resistance (Rm) plot
    ax4 = subplot(3,3,4);
    plot(wcp_corr_penn.Rm)
    title('Rm in MOhm')
    hold on; plot(movmean(wcp_corr_penn.Rm,7)); hold off
    % capacitace (Cm) plot
    ax5 = subplot(3,3,5);
    plot(wcp_corr_penn.Cm)
    title('Cm in pF')
    hold on; plot(movmean(wcp_corr_penn.Cm,7)); hold off
    % steady state current (Iss) plot
    ax6 = subplot(3,3,6);
    plot(wcp_corr_penn.Iss)
    title('Iss in pA')
    hold on; plot(movmean(wcp_corr_penn.Iss,7)); hold off
    % plot without test pulse
    ax7 = subplot(3,3,[7,9]);
    plot(x,YF,'Color',[0 0.4470 0.7410 0.2]) % to make lighter
    title('Recording w/ 1 kHz LPF')
    xlabel('Time (s)')
    ylabel('Current (pA)')
    hold on
    plot(x2,smooth(median_YF),'linewidth',2,'color',[0 0.4470 0.7410])
    hold off % not made lighter
    sgtitle('Compensated WCP Values by epoch')
    
    % set axis properties
    ax1.LineWidth = 2;
    ax1.Box = 'off';
    ax2.LineWidth = 2;
    ax2.Box = 'off';
    ax3.LineWidth = 2;
    ax3.Box = 'off';
    ax4.LineWidth = 2;
    ax4.Box = 'off';
    ax5.LineWidth = 2;
    ax5.Box = 'off';
    ax6.LineWidth = 2;
    ax6.Box = 'off';
    ax7.LineWidth = 2;
    ax7.Box = 'off';
    
    % save figure2
    a = split(file,'.');
    file_raw = append(char(a(1)),'_rscomp.pdf');
    saveas(gcf,file_raw)
    
    YF = [x' YF];
    figure; plot(YF(:,1), YF(:,2)); xlabel('time (s)'); ylabel('amp (pA)'); saveas(gcf,'tprm.fig');

    clear plots ax1 ax2 ax3 ax4 ax5 ax6 ax7 YF median_YF i x fraction ans file
    clear wcp_raw wcp_corr_penn
    clear S Param corr_p_splits x2 splits file_raw a
else
    disp("Plots disabled (to enable set plots = true)")
    clear('plots')
end

%% Define Functions


function y = rscomp (y, sampInt, fraction, R_s, C_m, V_hold, V_rev)

% Function for offline series resistance compensation

% y must a vector of the current trace (in pA)
% sampInt is sampling interval (in microseconds)
% fraction is the amount of compensation (1.0 = 100%)
% R_s is series resistance (in Mohm)
% C_m is cell capacitance (in pF)
% V_hold is holding potential (in mV)
% V_rev is reversal potential of the conductance (in mV)

numpoints = size(y,1);
numtraces = size(y,2);
if numtraces > 1
 error('This function only supports one trace/wave')
end
y = y * 1e-12;
sampInt = sampInt * 1e-6;
R_s = R_s * 1e6;
C_m = C_m * 1e-12;
tau = R_s * C_m;
V_hold = V_hold * 1e-3;
if nargin < 7
  V_rev = 0;
end
voltage = V_hold - V_rev;
tau_corr = R_s * C_m * (1 - fraction);

% First point: (we have to calculate this separately, because we need the value at i-1 below)
denominator = voltage - R_s * fraction * y(1);
if denominator ~= 0
  y(1) = y(1) * (voltage / denominator);
end

for i = 2:numpoints-1
  % this is the loop doing the correction for all other points
  % first calculate R_m for zero series resistance under the assumptions
  % that  U_m + U_Rs = const = voltage
  current = (y(i+1) + y(i)) / 2;  % The in between(mean) value
  derivative =  (y(i+1) - y(i)) / sampInt;
  denominator = current + tau * derivative;
  if denominator ~= 0
   R_m = (voltage - R_s * current) / denominator;  % calculate the true R_m
  else
    % Do nothing. Leave R_m as is
  end
  % Now calculate current for new series resitance
  denominator = (R_m + (1 - fraction) * R_s) * (1 + tau_corr / sampInt);
  if denominator ~= 0
   y(i) = tau_corr / (tau_corr + sampInt) * y(i-1) + voltage/denominator;
  else
   y(i) = y(i-1);  % old value
  end
end

y = y * 1e12;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function out = WCP(data,parameters)
%     input arguements:
%               data = x by n array of raw data where n is the number of
%                   sweeps or traces
%               paramaters = structure of parameters set in the experiment
%                   as defined at the start of this script
%     output arguemennts:
%               out = whole cell properties
%
%
% The baseline resting current, or holding current (Ih) is the trimmed
% average of the 10 seconds following the test pulse. Trimmed mean of 33%
% was used (Mendeleev, 1895)
Ih = zeros(size(data,2),1);
for i = 1:size(data,2)
    Ih(i) = trimmean(data(parameters.pulse_end:end,i),33,'floor');                   % in pA
end

% The series (or access) resistance is obtained my dividing the voltage step
% by the peak amplitude of the current transient (Ogden, 1994): Rs = V / Ip
base = zeros(size(data,2),1);
peak = zeros(size(data,2),1);
Ip = zeros(size(data,2),1);
Rs = zeros(size(data,2),1);
for i = 1:size(data,2)
    base(i) = trimmean(data(1:parameters.pulse_start,i),33,'floor');                              % baseline average (in pA)
    peak(i) = min(data(parameters.pulse_window,i));                                  % peak amplitude (negative)
    Ip(i) = base(i) -  peak(i);                                                % peak amplitude of the current transient (in pA)
    Rs(i) = ((parameters.voltage_step / 1000)/ Ip(i)*(10^12)) * 10^-6;           % in MOhms
end

% The input resistance is obtained by dividing the voltage step by the average amplitude of the steady-state current (Barbour, 2014): Rin = V / Iss
Iss = zeros(size(data,2),1);
Rin = zeros(size(data,2),1);
for i = 1:size(data,2)
    Iss(i) = trimmean(data(parameters.pulse_start:parameters.pulse_end,i),33,'floor') - base(i); % in pA
    Rin(i) = 0 - ((parameters.voltage_step / 1000) / Iss(i) *(10^12)) * 10^-6;    % in Mohms
end

% The cell membrane resistance is calculated by subtracting the series resistance from the input resistance (Barbour, 1994): Rm = Rin - Rs
Rm = zeros(size(data,2),1);
for i = 1:size(data,2)
    Rm(i) = Rin(i) - Rs(i);                                                 % in MOhms
end

% The cell membrane capacitance is ESTIMATED by dividing the transient charge by the size of the voltage-clamp step (Taylor et al. 2012): Cm = Q / V
Cm = zeros(size(data,2),1);
Q = zeros(size(data,2),1);
for i = 1:size(data,2)
    % calculate the estimated charge
    time = 0:5.0000e-05:0.01;                                              % in seconds
    Q(i) = trapz(time, ((base(i)-(data(parameters.pulse_window,i)))) * 1e-12); % charge in C
    % calculate the capacitance
    t = time(end)-time(1);
    Cm(i) = ((Q(i) - Iss(i) * t)   / (parameters.voltage_step / 1000));          % in pF
end

% assign variables to an output structure
out.Ih = Ih; % holding current in pA
out.Ip = Ip; % peak amplitude of the negative current transient in pA
out.base = base; % pre-pulse baseline in pA
out.peak = peak; % peak
out.Rs = Rs; % series resistance in Mohm
out.Iss = Iss; % intrapulse steady state current in pA
out.Rin = Rin; % input resistance in Mohm
out.Rm = Rm; % specific membrane resistance in Mohm
out.Cm = Cm; % capacitance in pF
out.Q = Q; % charge in C

% Other properties to be calculated later if required
%   The cell surface area is estimated by dividing the cell capacitance by thespecific cell capacitance, c (1.0 uF/cm^2; Gentet et al. 2000; Niebur, 2008):Area = Cm / c
%   The specific membrane resistance is calculated by multiplying the cell membrane resistance with the cell surface area: rho = Rm * Area
%   Users should be aware of the approximate nature of determining cell capacitance and derived parameters from the voltage-clamp step method (Golowasch, J. et al., 2009)

% References:
%    Barbour, B. (2014) Electronics for electrophysiologists. Microelectrode
%     Techniques workshop tutorial.
%     www.biologie.ens.fr/~barbour/electronics_for_electrophysiologists.pdf
%    Gentet, L.J., Stuart, G.J., and Clements, J.D. (2000) Direct measurement
%     of specific membrane capacitance in neurons. Biophys J. 79(1):314-320
%    Golowasch, J. et al. (2009) Membrane Capacitance Measurements Revisited:
%    Dependence of Capacitance Value on Measurement Method in Nonisopotential
%     Neurons. J Neurophysiol. 2009 Oct; 102(4): 2161-2175.
%    Niebur, E. (2008), Scholarpedia, 3(6):7166. doi:10.4249/scholarpedia.7166
%     www.scholarpedia.org/article/Electrical_properties_of_cell_membranes
%    (revision #13938, last accessed 30 April 2018)
%    Ogden, D. Chapter 16: Microelectrode electronics, in Ogden, D. (ed.)
%     Microelectrode Techniques. 1994. 2nd Edition. Cambridge: The Company
%     of Biologists Limited.
%    Taylor, A.L. (2012) What we talk about when we talk about capacitance
%     measured with the voltage-clamp step method J Comput Neurosci.
%     32(1):167-175

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [tprm,med,x,x2] = TestPulse_rm(data,parameters)

% remove test pulse and plot without
% calculate the pre pulse mean to fill the gaps
tprm_splits = data;
base = zeros(size(tprm_splits,2),1);
median_YF = zeros(size(tprm_splits,2),1);
for i = 1:size(tprm_splits,2)
    base(i) = trimmean(tprm_splits(1:parameters.pulse_start,i),33,'floor'); % baseline average (in pA)
    tprm_splits(parameters.pulse_start:(parameters.pulse_end+250),i) = base(i); % swap pulse for base
    median_YF(i) = median(tprm_splits(:,i)); % calculate the median of the splits
end
% concatenate tprm_splits
tprm_splits_conc = vertcat(tprm_splits(:));
% apply a filter to clean the look of the data
t = (0:size(tprm_splits_conc,1)-1)';
t = t./parameters.sample_rate;
YF = filter1(tprm_splits_conc, t, 0, 300);
x = linspace(0,(size(tprm_splits_conc,1)/20000),size(tprm_splits_conc,1));
%figure
%plot(x,YF,'Color',[0 0.4470 0.7410 0.2]) % to make lighter
%title('filter')
x2 = linspace(0,(size(tprm_splits_conc,1)/20000),size(median_YF,1));
%hold on
%plot(x2,smooth(median_YF),'linewidth',4,'color',[0 0.4470 0.7410]) % not made lighter
tprm = YF;
med = median_YF;
end