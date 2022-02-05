function out = WCP(data,parameters)
%     input arguements:
%               data = x by n array of raw data where n is the number of
%                   sweeps or traces (in pA)
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
    time = 0:1/parameters.sample_rate:0.01;                                              % in seconds
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
