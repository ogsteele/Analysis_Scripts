function NaACon = NaACon(Data)
for b = 1:37
    NaACon(b) = NaANormConductance(Data,b);
end
