% The cell membrane capacitance is ESTIMATED by dividing the transient charge by the size of the voltage-clamp step (Taylor et al. 2012): Cm = Q / V
Cm = zeros(size(splits,2),1);
Q = zeros(size(splits,2),1);
pre_pulse = zeros(size(splits,2),1);
post_Q_ind = zeros(size(splits,2),1);
pre_Q_ind = zeros(size(splits,2),1);
raw_upspike = zeros(1001,size(splits,2),1);
new_upspike = zeros(100001,size(splits,2),1);
for i = 1:size(splits,2)
    % calculate the pre pulse mean
    pre_pulse(i) = trimmean(splits(1:pulse_start,i),33,'floor');               % mean of the pre_pulse period
    % interpolate the trace to increase the data points
    % this increases the accuracy of the AUC for the charge
    raw_upspike(:,i) = splits(pulse_end:pulse_end+1000,i);                       % 1000 points = 20th of 1 second = 50 ms after the end of the test pulse
    raw_time = 0:1:length(raw_upspike(:,i))-1;                                   % create raw time array
    sample_fold = 100;                                                      % fold increase in sampling
    new_time = 0:1/sample_fold:length(raw_upspike(:,i))-1;                       % define new time arrayy
    new_upspike(:,i) = interp1(raw_time,raw_upspike(:,i),new_time);                   % 1D spline interpolation
    % find the location of the trace either side of the peak of the pulse
    % to allow us to calculate the AUC, and thus the charge
    [~,max_ind] = max(new_upspike(:,i));                                         % find the maximum point in the new trace
    post_Q_ind(i) = (find(new_upspike(max_ind:end,i) < pre_pulse(i),1)) + max_ind;              % pre pulse index and post pulse index
    pre_Q_ind(i) = find(new_upspike(1:max_ind,i) < pre_pulse(i),1,'last');          % find the last value lower than the pre_pulse average before the peak 
    % calculate the estimated charge
    Q(i) = trapz(new_upspike(pre_Q_ind(i):post_Q_ind(i),i)); % charge in coloumbs? What unit?
    % calculate the capacitance
    Cm(i) = Q(i) / (Param.voltage_step / 1000); % capacitance in Farads - what unit, need to confirm the above ... 
end