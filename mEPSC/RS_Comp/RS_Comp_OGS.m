%% Rs_Comp_OGS
% Script to load in ephys data, and compensate the series resistance
% changes observed during whole cell patch clamp recordings

% Note, change parameters as relevent

% Dependencies - ensure in path
% ephysIO (https://github.com/acp29/eventer/blob/master/base/ephysIO.m)

%% TO DO 
% legacy filetype in filename, split by delimiter
% optional low pass filter inclusion - not priority
% ensure continuity present with if statements
% format graphs too

%% Plots?
% toggle for plots
plots = false; % true / false (logical)
comp_penn = true; % true / false (logical)
comp_keine = false; % true / false (logical)
%% Parameters
% Amplifier / Sampling settings
Param.sample_rate = 20000; % in Hz
Param.amplifier_scale = 2e09; % actually 2e-09, however 2e09 reverses this
Param.amplifier_gain = 100;
Param.split_length = 10; % in seconds, frequency of test pulses
% Test pulse settings
Param.pulse_amp = -0.002; % in V
Param.Vh = -45; % holding voltage in mV
Param.commander_scale = 20; % in mV/V
Param.voltage_step = 2; % in mV
Param.pulse_duration = 10; % in ms
Param.pulse_points = (Param.pulse_duration*Param.sample_rate)/1000;               % convert from ms to data points
Param.pulse_start = 19100;
Param.pulse_end = Param.pulse_start + Param.pulse_points;                                     % note, this will only cover the downward transiet, not the upward transient
Param.pulse_window = Param.pulse_start:Param.pulse_end;
% Compensation settings
Param.Vrev = 0; % reversal potential in mV (close to zero for AMPAR/NMDAR)
Param.des_Rs = 10; % desired series resistance for recordings to be compensated to



%% Load in w/ ePhysIO 
% plus UI for script

% !! temporerily silenced whilst in development !!

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
    % save raw recording following the ephysIO guidance
    % To save a column-major XY array of electrophysiology data:
    %    ephysIO (filename,array,xunit,yunit,names,notes,datatype)
    % ephysIO((append('uncomp_',filename,'.phy')),S.array,S.xunit,S.yunit,S.names,S.notes,'int16')
    
end
% tidy workspace
clear('file', 'path','ans')


% For use during development - replacing the file dialog
%cd('/Users/ogsteele/Library/Mobile Documents/com~apple~CloudDocs/DPhil/Analysis/Analysis_Scripts/mEPSC')
%S = ephysIO('20201021_000.tdms');
%% Split the data into ten second waves
% convert the data in pA (like seen in nidaq_scope)
S.array(:,2) = S.array(:,2) * (Param.amplifier_scale * Param.amplifier_gain);
% calculate the length of the recording
length_raw = length(S.array(:,2));
length_seconds = length_raw/Param.sample_rate;
n_splits = length_seconds/Param.split_length;
% number of data points per split
length_raw_split = length_raw/n_splits;

% create list of start and end points determined by the sample rate, length
% of split, length of recording
% round down n_splits to not deal with incomplete sections and preallocate
start(:,1) = zeros(floor(n_splits),1);
finish(:,1) = zeros(floor(n_splits),1);
% define the starting numbers
start(1,1) = 1;
finish(1,1) = 1 + length_raw_split-1;
for i = 2:floor(n_splits)
    start(i,1) = finish(i-1,1) + 1;
    finish(i,1) = start(i,1) + length_raw_split-1;
end
    
% create a matrix of the split recording and preallocate
splits = zeros(length_raw_split,floor(n_splits));
for s = 1:floor(n_splits)
    splits(:,s) = S.array(start(s):finish(s),2);
end

% tidy up workspace
clear('start','finish','i','s','n_splits','length_seconds','length_raw', ...
    'length_raw_split')

%% Generate necessary whole cell paramaters
% locate the test pulse
wcp_raw = WCP(splits,Param); 

%% Compensation Penn
% adapted from rscomp_Penn.m
% 
if comp_penn == true
    disp("Comp_penn enabled (to disable set comp_penn = false)")
    % preallocate for speed 
    corr_p_splits = zeros(size(splits,1),size(splits,2));
    for i = 1:size(splits,2)
    corr_p_splits(:,i) = rscomp(splits,wcp_raw,Param,0.7,i); 
    end
else
    disp("Comp_penn disabled (to enable set comp_penn = true)")
end

% calculate whole cell properties from the compensated traces
wcp_corr_penn = WCP(corr_p_splits,Param);
%% Compensation Keine
% uses Christian Keines adaptation of Stephen Treynelis paper
% (https://www.sciencedirect.com/science/article/pii/S016502709800140X?via%3Dihub)

	% performs offline series resistance correction for recorded currents
	% Input:
	%
	% - data = recorded current trace, can be single vector or multiple recordings as array or cell
	% - Rs = uncompensated series resistance (Rs) during the recording in Ohm, i.e. if the Rs during the experiment was 10 MOhm and online compensated by 50% by
	% the amplifier, the remaining uncompensated Rs will be 5 MOhm (5e6 Ohm = 5 MOhm)
	% - Cm = Membrane capacitance during the recording in Farad (e.g. 10e-12 F = 10 pF)
	% - Vhold = holding potential during the recording in Volts (e.g. -0.06 V =  -60 mV)
	% - Vrev = reversal potential of the recorded current in Volts (e.g. 0.01V = 10 mV)
	% - SR = Sampling rate during the recordings (in Hz)
	% - [optional] fractionV = voltage error to be compensated [0-1] (e.g. 1 if voltage error should be fully compensated)
	% - [optional] fractionC = fraction capacitative filtering error to be compensated [0-1] (e.g. 1 capacitative filtering error should be fully compensated)
	% - [optional] fc = cutoff frequency for filter to smooth capacitative current correction (in Hz) (if omitted, fc is calculated from the sampling interval as fc = 1/(2 * pi * si))

	% Based on: "Traynelis SF (1998) Software-based correction of single
	% compartment series resistance errors. J Neurosci Methods 86:25–34."
	%
	% EXAMPLE: RsCorrection(data, Rs, Cm, Vhold, Vrev, SR, 'fractionV', 1, 'fractionC', 1, 'fc', 10e3)
    
    % adapted by OGS 03.13.20 to compensate to a set Rs Value by changing
    % the 'fractionV' value in the below code. 
if comp_keine == true
    % relative fraction V
    fractionV_calc = (Param.des_Rs ./ wcp_raw.Rs);
    fractionC_calc = fractionV_calc;
    disp("Comp_Keine enabled (to disable set comp_keine = false)")
    % run keine compensation in loop
    for x = 1:size(splits,2)
        Keine(x) = RsCorrection(splits(:,x), ...
                                (wcp_raw.Rs(x) * 1e6), ...
                                (wcp_raw.Cm(x) * 1e-12), ...
                                (Param.Vh * 1e-3), ...
                                (Param.Vrev * 1e-3), ...
                                Param.sample_rate, ...
                                'fractionV', fractionV_calc(x), ...
                                'fractionC', fractionC_calc(x), ...
                                'fc', 0);
    end
    
    % temporary diagnostic plots
    % arrange the compensated traces like splits
    corr_splits = zeros(length(Keine(1,1).dataRaw),size(Keine,2));
    for i = 1:size(Keine,2)
        corr_splits(:,i) = Keine(1,i).dataCorrected;
    end
    % calculate wcp for corrected traces
    wcp_corr_keine = WCP(corr_splits,Param); 
    % check for the peak amplitude of the upward test pulse spike
    figure
    plot(wcp_raw.Rs)
    hold on
    plot(wcp_corr_keine.Rs,'r')
    title('Blue = Rs raw, Orange = Rs corrected')
    % early trace corr
    figure
    plot(Keine(1,1).dataCorrected,'r')
    title('early corr')
    % late trace corr
    figure
    plot(Keine(1,size(Keine,2)).dataCorrected,'r')
    title('late corr')
    % early trace raw
    figure
    plot(Keine(1,1).dataRaw)
    title('early raw')
    % late trace raw
    figure
    plot(Keine(1,size(Keine,2)).dataRaw)
    title('late raw')
    % series resistance (Rs) plot
    figure
    plot(wcp_raw.Rs)
    title('Rs in MOhm')
    
else
    disp("Comp_keine disabled (to enable set comp_keine = true)")
end

% check for the peak amplitude of the upward test pulse spike

%% Plots
if plots == true
    disp("Plots enabled (to disable set plots = false)")
    % holding current (Ih) plot
    figure
    plot(Ih)
    title('Ih in pA')
    % series resistance (Rs) plot
    figure
    plot(Rs)
    title('Rs in MOhm')
    % raw trace plot 
    figure
    plot(S.array(:,2))
    title('Raw Trace')
    % input resistance (Rin) plot
    figure
    plot(Rin)
    title('Rin in MOhm')
    % memrane resistance (Rm) plot 
    figure
    plot(Rm)
    title('Rm in MOhm')
    % charge (Q) plot
        % figure
        % plot(Q)
        % title('Q in C')
    % capacitace (Cm) plot
    figure
    plot(Cm)
    title('Capacitance in pF')
    % steady state current (Iss) plot 
    figure
    plot(Iss)
    title('Iss in pA')
    
else
    disp("Plots disabled (to enable set plots = true)")
    clear('plots')
end 

%% Define Functions
function out = rscomp(data,properties,parameters,frac,epoch) 
%     input arguements:
%               data = x by n array of raw data
%               properties = structure of whole cell properties as
%                   generated by WCP()
%               paramaters = structure of parameters set in the experiment
%                   as defined at the start of this script
%               frac = compensation fraction (0-1) where 1 = 100%
%                   compensation
%               epoch = trace selection
%     output arguemennts:
%               out = compensated data
%
%
%
%     This script is my Matlab/Octave adaptation of Erwin Neher's Igor functions:
%       RunSeriesResComp
%       SeriesresistanceComp
%     These were freely available in in Proc02_Apr4Web.ipf from Erwin Neher's webpage:
%      http://www3.mpibpc.mpg.de/groups/neher/index.php?page=software
%      (last accessed: 01 July 2014)
%
%     The function replaces current traces by their series-compensated version;
%     the value at i is replaced by the average at i and i+1
%     R_s is in ohms, C_m in Farads, fraction is the fraction to be compensated
%     if R_s was 5 MOhm in the experiment and if it was 50% hardware compensated
%     then R_s = 2.5e6 has to be entered and f=1 for a complete overall compensation
%     The routine, similarly to that published by Traynelis J. Neurosc. Meth. 86:25,
%     compensates the frequency response, assuming a single R_s*C_m - time constant
%     (at constant V-hold) and a three component equivalent circuit for the pipette cell
%     assembly with  R_s, C_m, R_m
%
%     Theory: V_h = R_s*I+U_m;
%     I_r is membrane resistive current,
%     I is total current
%     I_r = I-I_c = I-C_m*dU_m/dt = I+C_m*R_s*dI/dt (because R_s*I+U_m = const.)
%     G_m=I_r/ U_m = (I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
%     For complete correction (fraction = 1) : I_corr = V_h*(I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
%
%     rscomp v1.0 (last updated: 01/06/2014)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/

% legacy dependancies 
R_s = properties.Rs(epoch) * 1e6;                                           % series resistance in Ohms
C_m = properties.Cm(epoch) * 1e-12;                                         % capacitance in Farads
tau = R_s * C_m;                                                            % time constant
V_hold = parameters.Vh * 1e3;                                               % holding potential in V
V_reversal = parameters.Vrev * 1e3;                                         % reversal potential in V
voltage = V_hold - V_reversal;                                              % voltage difference between the holding and reversal in V
fraction = frac;                                                            % compensation fraction
tau_corr = R_s * C_m * (1 - fraction);                                      % tau correction
y = data(:,epoch);                                                          % data from input arguement
numpoints = length(data);                                                   % number of data points
sampInt = 1/parameters.sample_rate;                                         % sampling interval in s

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
 out = y; % compensated data?
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