function Out = it_wrapper(filename)

Location = readtable(filename);

    for i = 1:height(Location)
        Out(i) = importtest(char(Location{i,1}));
    end
    
end