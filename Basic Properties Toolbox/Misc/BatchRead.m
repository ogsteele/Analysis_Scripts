function [a,b,c,d] = BatchRead

abf = dir('*.abf'); 
numfiles = length(abf);
mydata = cell(1, numfiles);

for k = 1:numfiles 
  mydata{k} = abfload(abf(k).name); 
end
[a] = mydata{1,1};
[b] = mydata{1,2};
[c] = mydata{1,3};
[d] = mydata{1,4};

end
 