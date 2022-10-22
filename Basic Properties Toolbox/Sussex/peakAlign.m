function alignedArray = peakAlign(y,time)
%
% input arguments
%   array = x * y array output from ephysIO
%   time = [a,b] of time points in data points
%
% output arguments
%   alignedArray = peak aligned array of peak within time
%
% example usage
%   S = ephysIO(clampfile)
%   x = S.array(:,1)
%   y = S.array(:,2:end)
%   time = [1000,4000]
%   alignedArray = peakAlign(y,time)
%
% Note
%   Data loss is quite aggressive, but works for this usage. Consider
%   improving in future if use is to become more ubiquitous in code.
%% CODE

% locate minimum indexes within region of interest
 [~,ind] = min(y(time(1):time(2),:));
corrInd = ind + time(1); % correct for drift associated with peak detection

% determine the peak differences relative to the smallest
[val,~] = min(corrInd); % earliest peak and wave number
indDiff = (corrInd - val)+1; % differences in each peak index

% trim the start of the arrays, and concatenate arrays
alignedArray = zeros(sum(time)*2,size(y,2));
for i = 1:size(y,2)
    alignedArray(:,i) = y(indDiff(i):(indDiff(i)+(sum(time)*2))-1,i);
end

end
