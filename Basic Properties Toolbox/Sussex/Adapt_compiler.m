


for i = 1:size(compiled,2)
    filepath = [compiled(i).filepath,'Clamp1.ma'];
    adaptation(i) = SpikeAdapt(300,1,filepath);
end


Genotype = {compiled.genotype};

Ketamine = {compiled.ketamine};
Memantine = {compiled.memantine};

for c = 1:size(Genotype,2)
    if Ketamine(c) == "No" && Memantine(c) == "No"
        condition(c) = 'Control';
    elseif Ketamine(c) == "Yes" && Memantine(c) == "No"
        condition(c) = 'Ketamine';
    elseif Ketamine(c) == "No" && Memantine(c) == "Yes"
        condition(c) = 'Memantine';
    end
end

Genotype = Genotype';
Genotype = cell2mat(Genotype);
Condition = condition';

for km = 1:size(adaptation,2)
    adaptation(km).genotype = Genotype(km,:);
    adaptation(km).condition = char(Condition(km));
    adaptation(km).numSpikes = compiled(km).numSpikes(26);
    if adaptation(km).numSpikes > 3
        adaptation(km).exclude = 0;
    else 
        adaptation(km).exclude = 1;
    end
end

%% Split into condition
E4_control = adaptation(strcmp({adaptation.genotype}, 'APOE4') & strcmp({adaptation.condition}, 'Control'))
for i=1:numel(E4_control)
   if ~isequal(logfile(i).exclude,'1')
      E4_control(i)=[];
   end
end
E3_control = adaptation(strcmp({adaptation.genotype}, 'APOE3') & strcmp({adaptation.condition}, 'Control'))
E4_ketamine = adaptation(strcmp({adaptation.genotype}, 'APOE4') & strcmp({adaptation.condition}, 'Ketamine'))
E3_ketamine = adaptation(strcmp({adaptation.genotype}, 'APOE3') & strcmp({adaptation.condition}, 'Ketamine'))
E4_memantine = adaptation(strcmp({adaptation.genotype}, 'APOE4') & strcmp({adaptation.condition}, 'Memantine'))
E3_memantine = adaptation(strcmp({adaptation.genotype}, 'APOE3') & strcmp({adaptation.condition}, 'Memantine'))

% manually delete the crap ones now

E3C_Hz = NaN(size(E3_control,2),20);
for row = 1:size(E3_control,2)
    for collumn = 1:size(E3_control(row).peakHz,1)
    E3C_Hz(row,collumn) = E3_control(row).peakHz(collumn);
    end
end

E4C_Hz = NaN(size(E4_control,2),20);
for row = 1:size(E4_control,2)
    for collumn = 1:size(E4_control(row).peakHz,1)
    E4C_Hz(row,collumn) = E4_control(row).peakHz(collumn);
    end
end

E4K_Hz = NaN(size(E4_ketamine,2),20);
for row = 1:size(E4_ketamine,2)
    for collumn = 1:size(E4_ketamine(row).peakHz,1)
    E4K_Hz(row,collumn) = E4_ketamine(row).peakHz(collumn);
    end
end

E3K_Hz = NaN(size(E3_ketamine,2),20);
for row = 1:size(E3_ketamine,2)
    for collumn = 1:size(E3_ketamine(row).peakHz,1)
    E3K_Hz(row,collumn) = E3_ketamine(row).peakHz(collumn);
    end
end

E3M_Hz = NaN(size(E3_memantine,2),20);
for row = 1:size(E3_memantine,2)
    for collumn = 1:size(E3_memantine(row).peakHz,1)
    E3M_Hz(row,collumn) = E3_memantine(row).peakHz(collumn);
    end
end

E4M_Hz = NaN(size(E4_memantine,2),20);
for row = 1:size(E4_memantine,2)
    for collumn = 1:size(E4_memantine(row).peakHz,1)
    E4M_Hz(row,collumn) = E4_memantine(row).peakHz(collumn);
    end
end
