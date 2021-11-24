% cd = change directory
% char = converts cell to character array
% {} = extracts singular data from a table
% nnz = number of non zeroes in a matrix, useful for counting
% sum = counts the total. 
% mean includes the zeros, therefore sum(a)/nnz(a) returns accurate mean
% std = standard deviation
% ([Data(:).x]) returns all values in the field, x


% MeanX = mean([Data(:).MinNaADensity], 'omitnan')
% nX = sum(not(isnan([Data(:).MinNaADensity])),2)
% SEMX = (nanstd([Data(:).X])/(sqrt(nX))
% StDevX = nanstd([Data(:).X])
