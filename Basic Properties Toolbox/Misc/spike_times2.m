function [N ,out1] = spike_times2(trace,threshold1)
%   This function detects and locates the time points of action potentials in a trace of 
%   membrane potential as a function of time in a neuron. The trace should represent
%   a current clamp recording from a neuron.
%   Input: 
%   "trace" is the membrane voltage array of the neuron
%   "Theshold" is the value for the spike to cross to be detected.
%   Output:
%   N - the number of spikes detected
%   out1 - array of sample numbers (i.e. times) of spikes detected
%
%   Example:
% 
%   [numberOfSpikes, spikeTimes] = spike_times2(d,-10);
%
%   assuming "d" contains trace loaded previously, and -10 mV was chosen as
%   action potential threshold. 
% 
%   Original script by: 
%   Rune W. Berg 2006 / rune@berg-lab.net / www.berg-lab.net
%   Modified by: Rune W. Berg, May 2015
%   Modified by: Roman E. Popov, June 2016. (rpopov@uvm.edu)

 gim=trace;

    clear('set_crossgi')
    % find indices of samples above the threshold:
    set_crossgi=find(gim(1:end) > threshold1)  ; 
    clear('index_shift_neggi');clear('index_shift_pos');
    
% if there is at least one sample above the threshold:
if isempty(set_crossgi) < 1  
    clear('set_cross_plusgi');clear('set_cross_minus')
    % first index posgi = 1st sample above the threshold:
    index_shift_posgi(1)=min(set_crossgi);
    % set last index neggi to the last sample above the threshold:
    index_shift_neggi(length(set_crossgi))=max(set_crossgi);

    for i=1:length(set_crossgi)-1
        % this line detects when there is a discontinuous jump in the
        % indices, thus -> shift from one spike to another. Example: 
        % a) continuous indices: set_crossgi(i+1)=721 > set_crossgi(i)+1 =
        % 720 + 1; therefore, 721 = 721, not true, 0, skip the clause.
        % b) discontinuous, set_crossgi(i+1)= 1000 > set_crossgi(i)+1 =
        % 721+1 = 1000 > 722, true -> record the indices
        if set_crossgi(i+1) > set_crossgi(i)+1 ; 
            % index of the right end of the spike interval
            index_shift_posgi(i+1)=i;
            index_shift_neggi(i)=i;
        end
    end
    % ^^^ Code extracts indices of the "above treshold" samples

    %Identifying up and down slopes:
    % Here indices of the truncated "above the threshold" trace are
    % converted into sample indices of the original trace:
    set_cross_plusgi=  set_crossgi(find(index_shift_posgi));   % find(x) returns nonzero arguments.
    set_cross_minusgi=  set_crossgi(find(index_shift_neggi));   % find(x) returns nonzero arguments.
    set_cross_minusgi(length(set_cross_plusgi))= set_crossgi(end);
    

    nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows.
    

    % Getting the spike coords
    spikemax = zeros(1, nspikes);
    spikemaxCorrected = zeros(1, nspikes);
    
    for i=1: nspikes
            % identify a potential spike:
            spikemax(i)=min(find(gim(set_cross_plusgi(i):set_cross_minusgi(i)) == max(gim(set_cross_plusgi(i):set_cross_minusgi(i))))) +set_cross_plusgi(i)-1;
            % from second spike and on, check that it is not too close or
            % even within the prev. spike inteval:
            if i>1
                % to get the average, set the corrected:
                spikemaxCorrected(i) = spikemax(i);

                % avg. time between spikes:
                avgIntraSpikeTime = mean(diff(spikemaxCorrected(spikemaxCorrected>0)));
                % check if the spike detected is within the previous
                % spike's interval:
                if ismember(spikemax(i),set_cross_plusgi(i-1):set_cross_minusgi(i-1))
                      % decide which spike to keep (based on trace values -> higher value means true peak)
                      if trace(spikemax(i-1))>trace(spikemax(i))
                          % previous spike was higher, discarding the
                          % current value:
                          spikemaxCorrected(i) = 0;
                      else
                          % current spike is higher, discarding the
                          % previous and keeping the current one: 
                          spikemaxCorrected(i-1) = 0;
                          spikemaxCorrected(i) = spikemax(i);
                      end
                % check if current spike is too close to the previous to be a real one.
                % True spike should be at least later than 50% of average intra-spike interval      
                elseif (spikemax(i) - spikemax(i-1)) < 0.5*avgIntraSpikeTime
                      % Similiar to above, keep the spike with higher amplitude:   
                      if trace(spikemax(i-1))>trace(spikemax(i))
                          spikemaxCorrected(i) = 0;
                      else
                          spikemaxCorrected(i-1) = 0;
                          spikemaxCorrected(i) = spikemax(i);
                      end
%                 else
%                     spikemaxCorrected(i) = spikemax(i);                    
                end
            else
                % Assuming that the 1st spike is always correct 
                % (this will be verified below):
                spikemaxCorrected(i) = spikemax(i);
            end
    end
else
    spikemax=[];
    spikemaxCorrected = [];
    display('no spikes in trace')
end

% remove erroneous spikes:
clear spikemax
spikemax = spikemaxCorrected(spikemaxCorrected>0);
% run one more check (this will remove erroneous spikes when there were 
% several erroneous spikes in a row and if the first one was incorrent):
idx2remove = [];
for peak = 2:length(spikemax)
    diffBtwSpikes = spikemax(peak) - spikemax(peak-1);
    if diffBtwSpikes < 0.5 * avgIntraSpikeTime
        if trace(spikemax(peak)) < trace(spikemax(peak-1))
            idx2remove = [idx2remove peak];
        else
            idx2remove = [idx2remove peak-1];
        end
    end
end
spikemax(idx2remove)=[];
% plot:
figure; plot(trace); hold on; plot(spikemax, trace(spikemax),'or');hold off

 
N=length(spikemax) ;

out1=spikemax;


