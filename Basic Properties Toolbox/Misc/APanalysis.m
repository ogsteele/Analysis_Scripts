function APanalysis(filename,threshold)
% APanalysis a quick tool to determine the number of spikes, their
% latencies and amplitudes, spiking frequency, action potential threshold,
% 50% widths (+errors in this estimate due to sampling). 
% Thresholds are hand-picked, which takes much longer than automatic, but
% is much more reliable.
%
% Argument list is:
% filename - full name of the file with .abf in the end: '123456.abf'
% threshold - minimal voltage for a spike to be considered as action
% potential, not just membrane potential fluctuation. Set this based on
% your recordings (e.g., if all of you action potentials are at least 
% above 5 mV, then pass threshold = 5). By default, threshold = 0.
%
% Manual input will wait for the same number of inputs as the number of
% spikes identified in the current trace. "Enter" key should terminate
% manual input earlier, but might cause error in this script. This issue
% may be addressed later.
% 
% To close all the plots when done: 
% >> close all
% 
% Script version 1.2, updated June 2016. 
% Created and modified by Roman E. Popov (rpopov@uvm.edu)
 

if nargin < 2
    height = 0;
else
    height = threshold;
end

% Load data from .abf file. Using newer version abf2load:
[d,si] = abfload(filename,'start',1,'stop','e');

totalSweeps = size(d,3);
% Initialize the var for frequency of each sweep:
freq = zeros(1,totalSweeps);  

for sweep=1:totalSweeps
    %%SPIKE DETECTION
    [count,peak_idx] = spike_times2(d(:,1,sweep),height);  
    
    %Frequency calculation
    if isnan(count)
        freq(sweep) = NaN;
    else
        secs_recorded = size(d,1)*si/10^6; %total # of samples * us in on 
        % sample / 1,000,000 - since the microseconds are converted to seconds for frequency calculation
        freq(sweep) = count/secs_recorded; 
    end
    
    % Initialize the var for amplitudes and amplitude times:
    amps = zeros(count,2);% the first value is the amplitude of the peak in uV, the second - the time of the amplitude in sec
    % Header for the .xls export file:
    header=[{'Peak Amplitudes, mV'} {'Amplitude Time, msec'} {'AP Freq, HZ'} ...
        {'Threshold, mV'} {'Threshold reached at, msec'} {'AP 50% width, msec'}...
        {'Error in AP 50% width (LEFT), %'} {'Error in AP 50% width (RIGHT), %'}];
    for l=1:count
        % The amplitudes are calculated later, relative to
        % threshold voltage
        amps(l,2) = peak_idx(l)*si*1e-3;   %Amplitude times converted to milliseconds
    end
    % if there were no spikes detected:
    if isempty(amps)
        amps(:,1)=NaN;
        amps(:,2)=NaN;
    end
    
    %%Work through each spike: 
    if ~isempty(peak_idx)
        side_length = floor(size(d,1)/(2*count));%Get the length of the side of a peak
        % initialize flags for begin/end that are shorter than the
        % side_length:
        short_begin = 0;
        short_end = 0;
        
        % Check if there is a long period without AP in the begin/end of
        % the trace:
        trace_begin = size(d(1:peak_idx(1),1,sweep),1);
        trace_end = size(d(peak_idx(end):end,1,sweep),1);
        
        if trace_begin > side_length*2  
        % if there is a silent period in 
        % the beginning of the trace that is longer than twice the half of 
        % AP width then it means that the side_length was calculated
        % improperly, and silent period should be truncated for proper
        % estimation of side_length:
            begin_shift = trace_begin - side_length;
            trunc_d = d(1+begin_shift:end,1,sweep);
            % recalculate side_length of an AP:
            side_length = floor(size(trunc_d,1)/(2*count));
        elseif trace_begin < side_length
            short_begin = 1;
        end
        
        % Check the same for the end of the trace:
        if trace_end > side_length*2 
            trunc_d = d(1:end - side_length,1,sweep);
            % recalculate side_length of an AP:
            side_length = floor(size(trunc_d,1)/(2*count));
        elseif trace_end < side_length
            short_end = 1;
        end
      
        location1 = zeros(1, count);%For width determination, left side
        location2 = zeros(1, count);%For width determination, right side
        
        APwidth = zeros(count,1);% AP width at 50% of amplitude
        x = zeros(1,count); % threshold time
        y = zeros(1,count); % threshold voltage
        
        % Due to sampling, it is impossible to find the time step at which
        % voltage is exactly equal 50% value. Thus, next one to it is
        % chosen. To assess how much this estimate is off relative to the
        % amplitude itself, these parameters are calculated:
        left_error = zeros(1,count); 
        right_error = zeros(1,count);

        % Calculate the width @50%
        prev_min_idx = 0;
        for k=1:count
            % make figure fullscreen to ease picking the thresholds:
            screenSize = get(0,'Screensize');
            set(gcf, 'Position', screenSize); % Maximize figure
            % get the point from the figure:
            % [x(k),y(k)] = ginput(1);
            [x_temp,y_temp] = ginput(1);
            if ~isempty(x_temp)&&~isempty(y_temp)
                x(k) = x_temp;
                y(k) = y_temp;
                % plot the chosen point in the figure itself (to confirm it was picked)
                hold on
                plot(x(k),y(k),'ko')
            else
                x(k) = NaN;
                y(k) = NaN;
            end

                    
            % Width determination
            % On each of the slopes I find the last point that is less 50% (left 
            % side), and first point <50% (right side).

            if k == 1
                left_border = 1;       
            else                            
                left_border = prev_min_idx;% this value is from the previous cycle
            end

            if short_end
                if k == count
                    right_border = size(d,1); % number of elements = index of the last element; to avoid going outside the array with index too large
                end
            else
                [~,right_border] = min(d(peak_idx(k):peak_idx(k)+side_length,1,sweep));
                right_border = right_border + peak_idx(k);
            end
            % The afterhyperpolarization value and location is used to
            % separate action potentials:
            [~,min_idx] = min(d(peak_idx(k):right_border,1,sweep));
            max_val = max(d(peak_idx(k):right_border,1,sweep));
            perc50val = y(k) + ((max_val-y(k))/2); % take 50% value and add it to this peak's minimum
            location1(k) = find(d(left_border:peak_idx(k),1,sweep) < perc50val,1,'last');%Since the AP is rising here I am looking for the last point before voltage reaches 50%
            location2(k) = find(d(peak_idx(k):right_border,1,sweep) < perc50val,1,'first');%Conversely, I am looking for the first sign that voltage dropped below 50% here, on the right side of the slope of AP

            
            % shift into absolute trace coordinates:
            location1(k) = location1(k) + left_border-1;
            location2(k) = peak_idx(k) + location2(k);
            
            % Check if the neighboring sample is closer to 50% voltage:
            if abs(perc50val - d(location1(k),1,sweep)) > abs(perc50val - d(location1(k)+1,1,sweep))
                location1(k) = location1(k) + 1;
            end
            
            if abs(perc50val - d(location2(k),1,sweep)) > abs(perc50val - d(location2(k)-1,1,sweep))
                location2(k) = location2(k) - 1;
            end
            
            % Assess how much off in voltage are the actual samples chosen
            % for AP width calculation. Negative values mean the estimate
            % is below the real value, positive - above:
            left_offset = -(perc50val - d(location1(k),1,sweep));
            left_error(k) = 100*left_offset/(max_val-y(k));

            right_offset = -(perc50val - d(location2(k),1,sweep));
            right_error(k) = 100*right_offset/(max_val-y(k));
  
            APwidth(k) = (location2(k)-location1(k))*si*1e-3;%Convert the width into milliseconds here as well
            amps(k,1) = max_val - y(k);%Amplitudes

            % Place the points used for APwidth estimation:
            % hold on    
            plot(location1(k),d(location1(k),1,sweep),'co');
            plot(location2(k),d(location2(k),1,sweep),'go');
            % add a 50% line for each spike
            hold off

            legend('Trace','AP Peaks','Picked Threshold','Left 50%', 'Right 50%')
            % Store the index of the lowest point of the current spike 
            % (usually, it is the point of afterhyperpolarization) to start
            % off analysis of the next spike:
            prev_min_idx = min_idx;

        end
   
        %Export of threshold voltages and AP halfwidths
        %header=[{'Peak Amplitudes, mV'} {'Peak Times, msec'} {'Freq, HZ'} {'Threshold, mV'} {'AP onset, msec'} {'AP 50% width, msec'}];
    xlswrite([filename '.xls'],header,['Sweep ' int2str(sweep)],'A1');    
    xlswrite([filename '.xls'],amps(:,1),['Sweep ' int2str(sweep)],'A2');
    xlswrite([filename '.xls'],amps(:,2),['Sweep ' int2str(sweep)],'B2');
    xlswrite([filename '.xls'],freq(sweep),['Sweep ' int2str(sweep)],'C2');
    xlswrite([filename '.xls'],y',['Sweep ' int2str(sweep)],'D2');
    xlswrite([filename '.xls'],(x*si*1e-3)',['Sweep ' int2str(sweep)],'E2');
    xlswrite([filename '.xls'],APwidth,['Sweep ' int2str(sweep)],'F2');
    xlswrite([filename '.xls'],left_error',['Sweep ' int2str(sweep)],'G2');
    xlswrite([filename '.xls'],right_error',['Sweep ' int2str(sweep)],'H2');
    
     end
    
end
% Possibly useful graphs of averaged statistics for the whole recording

% figure;
% plot(freq,'rx');
% title('Frequency of AP Firing in Every Sweep');
% 
% figure;
% plot(avg_thresh_voltage,'rx');
% title('Average Threshold Voltage in Every Sweep');
% 
% figure;
% plot(avg_width,'rx');
% title('Average AP width @ 50% in Every Sweep');


