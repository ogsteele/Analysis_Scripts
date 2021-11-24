function Density(Data)
for i = 1:length(Data)
    MeanNaAConductance(i) = Data(i).NaANormConductance(1);
end
mean(MeanNaAConductance, 'omitnan')
end