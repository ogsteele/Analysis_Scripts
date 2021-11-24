function [C_out,ia_out,ic_out] = testest(csvfilename)

[ChannelID,ChannelLabel,WellID] = import_csv_mea(csvfilename);

[C,ia,ic] = unique(WellID,'stable');

for x = 1:((length(C))-1)
        ia(x,2) = ia((x+1),1) - 1;
end
    ia(length(C),2) = length(ic);
    
C_out = C;
ia_out = ia;
ic_out = ic;
end