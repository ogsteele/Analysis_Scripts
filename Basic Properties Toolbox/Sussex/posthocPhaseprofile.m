function [mV,dvdt] = posthocPhaseprofile(condition)
% https://www.nature.com/articles/nrn2148
% https://www.cell.com/biophysj/pdf/S0006-3495(63)86827-7.pdf
% Tempted not to be too distracted by the phase plane plots 
% however as Bean et al suggest they're less informative when
% recorded from cell bodies due the confounding presence of the
% dendritic tree / axon. Not to mention these are also action 
% potentials after current injection, which further 
% complicates it


% post hoc phase profile of the AP script
% to be included in the IStep_compiler graph

% set the detection zones
S.xdiff = condition(1).ephysIO.xdiff;
detection = [2*1e4,6*1e4];
detStart = round(detection(1)/S.xdiff);

% set the winSize_ms and LPF_Hz variables
winSize_ms = 25;
winSize_ms = (winSize_ms*1e-3)/S.xdiff; % (ms to s)to data points
LPF_Hz = 330;

% figure; % figure when debugging
% construct for loop to run process
for i = 1:size(condition,2)
    % determine the wavenumber of first AP
    locs = condition(i).locs;
        C = locs;
            idx = ~cellfun('isempty',C);
            outs = zeros(size(C));
            outs(idx) = cellfun(@(v)v(1),C(idx));
            for lp = 1:size(outs,1)
                if outs(lp) > 0
                outs(lp) = outs(lp) + detStart; % account for the detection zone
                end
            end
            logicalIndexes =  outs < 27500 & outs > 1;
    wavenum_first(i) = (size(locs,2)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
   
    % detemine AP window for wave of interest
    warning('off','MATLAB:colon:nonIntegerIndex')
    mV(:,i) = condition(i).waves(outs(wavenum_first(i))-(winSize_ms/2):outs(wavenum_first(i))+(winSize_ms/2),wavenum_first(i))*1000;
    % calcuate dv/dt
    dvdt(:,i) = gradient(mV(:,i));
    warning('on','MATLAB:colon:nonIntegerIndex')
    %hold on; plot(mV(:,i)); legend() % plot when debugging
end
