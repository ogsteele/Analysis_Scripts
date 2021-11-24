function NaANormConductance(Data)

for i = 1:length(Data)
    MeanNaAConductance(i) = Data(i).NaANormConductance(b);
end
mean(MeanNaAConductance, 'omitnan')
end