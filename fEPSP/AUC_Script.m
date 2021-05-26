% clear workspace and close figures
clear
close

% cd to wd and load in raw data csv
cd('/Users/ogsteele/OneDrive - University of Sussex/R/Recombinant Washon');
data_Raw = csvread('Base.csv',5,1); 

% preallocate for speed
AUC = zeros(1,size(data_Raw,2));
base = zeros(1,size(data_Raw,2));

% for loop for the AUC calculation
for i = 1:size(data_Raw,2)
    % calculate the baseline as the mean of the first 20 minutes and populate
    % as a line
    base(i) = mean(data_Raw(1:40,(i)));
    baseline = zeros(80,1);
    baseline(1:end) = base(i);

    % calculate the AUC of trace relative to baseline AUC
    AUC_total = trapz(data_Raw(41:end-1,(i)));
    AUC_baseline = trapz(baseline);
    AUC(i) = diff([AUC_baseline, AUC_total]);
end

save('AUC','AUC','-ascii','-double','-tabs')
% Prism_AUC = [9.33887271124997	3.91749407549997	10.7350586537500 ...
%     -6.11522550027500	9.08176787562497		-6.11522541750001 ...
%     32.5752971452500	25.5661550042500	8.82574888525002	...
%     86.1681668917501	32.5083694597500	5.87538577374995	...
%     6.22785920549995	2.70889998475003		25.4802446990000...
%     3.52637093525000	3.40665558350003	-2.99721720349995	...
%     38.0230935935000	9.08176784300002	33.2177961210000	...
%     -2.07046088150000	0.749082624499994		2.88366052400000...
%     0.818443234749994	-2.41445407349997	20.1191484365000	...
%     8.36529493875000	42.7154779000000	30.7538185687500	...
%     54.6241828487500	35.6307537237500	22.4826775255000	...
%     -4.94776330424998	2.79868170275002	47.4061559380000	...
%     -3.38768266775000	18.5872258080000	41.4591012730000	...
%     5.35460350625000	13.6821841390000	14.5717873990000	...
%     15.3404476020000];
% 
% % calculate the difference between AUC calculations
% Difference = Prism_AUC - AUC;
% Difference;
% avg_diff = mean(Difference)
% [h,p,ci,stats] = ttest(Prism_AUC,AUC)