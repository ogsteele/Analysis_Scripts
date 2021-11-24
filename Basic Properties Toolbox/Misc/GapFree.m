function [Vm] = GapFree(filenameabf)
%GapFree 
%   abfload loads gapfree recording in abf format
gapfree = abfload(filenameabf);
%   mean of the first 10ms to provide Vm
Vm = mean(gapfree(1:200000));
end
