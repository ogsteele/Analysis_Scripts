function y = rscomp (y, sampInt, fraction, R_s, C_m, V_hold, V_rev)

% Function for offline series resistance compensation

% y must a vector of the current trace (in pA)
% sampInt is sampling interval (in microseconds)
% fraction is the amount of compensation (1.0 = 100%)
% R_s is series resistance (in Mohm)
% C_m is cell time constant C_m (in ms)
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
