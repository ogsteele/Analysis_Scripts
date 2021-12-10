[file,path] = uigetfile;
S = ephysIO(fullfile(path,file));
% preallocate pks, locs, w, p
warning('off','signal:findpeaks:largeMinPeakHeight');
for i = 2:size(S.array,2)
    [pks{i-1},locs{i-1},w{i-1},p{i-1}] = findpeaks(S.array(:,i),...
        'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',0.02*10^4);
    numSpikes(i-1) = size(pks{i-1},1);
end
% determine rheobase as the first amount of current to induce APs in the
% first 25% of the current step rather than the first current step value to
% elcit any AP at all. 