
[August02_WT_DIV10,August02_KO_DIV10,August02_MIX_DIV10] = SplitTestTwo('12-08-18 DIV10 G30-KO-MIX -1.csv','12-08-18 DIV10 G30-KO-MIX-2.csv');
[August02_WT_DIV12,August02_KO_DIV12,August02_MIX_DIV12] = SplitTestTwo('14-08-18 DIV12 G30-KO-MIX -1.csv','14-08-18 DIV12 G30-KO-MIX -2.csv');
[August02_WT_DIV14,August02_KO_DIV14,August02_MIX_DIV14] = SplitTestTwo('16-08-18 DIV14 G30-KO-MIX -1.csv','16-08-18 DIV14 G30-KO-MIX -2.csv');
[August02_WT_DIV16,August02_KO_DIV16,August02_MIX_DIV16] = SplitTestTwo('18-08-18 DIV16 G30-KO-MIX -1.csv','18-08-18 DIV16 G30-KO-MIX -1.csv');

% Concatenate structures
August02_WT = [August02_WT_DIV10;August02_WT_DIV12;August02_WT_DIV14;August02_WT_DIV16;];
August02_KO = [August02_KO_DIV10;August02_KO_DIV12;August02_KO_DIV14;August02_KO_DIV16;];
August02_MIX = [August02_MIX_DIV10;August02_MIX_DIV12;August02_MIX_DIV14;August02_MIX_DIV16;];