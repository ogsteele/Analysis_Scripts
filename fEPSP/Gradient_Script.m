% Script to calculate the slope of the raw data before application
% and after application of ApoE

% clear workspace and close figures
clear
close

% cd to wd and load in raw data csv
cd('/Users/ogsteele/OneDrive - University of Sussex/R/Recombinant Washon');
data_Raw = csvread('Base.csv',5,1); 

% preallocate for speed
p = zeros(43,2);
q = zeros(43,2);

% calculate the coefficients of linear interpolation - baseline
base_values = data_Raw(1:40,:);
for i = 1:size(data_Raw,2)    
    p(i,:) = polyfit((1:40)',base_values(:,i),1);
    Gradient_Base = (p(:,1)) * 2;
end

% calculate the coefficients of linear interpolation - post application
post_values = data_Raw(41:end-1,:);
for n = 1:size(data_Raw,2)    
    q(n,:) = polyfit((1:80)',post_values(:,n),1);
    Gradient_Post = (q(:,1)) * 2;
end

% calculate ratio of early to late gradients 
Gradient_Ratio = Gradient_Post ./ Gradient_Base;

% calculate percent difference 
Gradient_Change = (Gradient_Ratio - 1)*100; 

% calculate difference of early to late gradients
Gradient_Diff = Gradient_Post - Gradient_Base;

save('Gradient_Ratios','Gradient_Ratio','-ascii','-double','-tabs')
save('Gradient_Change','Gradient_Change','-ascii','-double','-tabs')
%% Stats for early

% % Slope values taken from Prism
% Prism_Gradient_Early = [0.001890	0.01508	-0.009108	-0.003660	0.007094 ...
%     -0.003660	0.02254	0.03167	0.03837	0.01535	0.01901	-0.01145	...
%     0.006074	0.01065		-0.01230	0.001760	0.002126	...
%     0.04561	0.004576	0.007094	0.02199	0.003430	-0.01299...
%     -0.0007358	0.01441	0.03150	0.01058	0.01853	0.04548	0.005008...
%     0.03119	0.02573	0.01014		0.009841	0.01932	0.02310	-0.001633...
%     0.03063	0.01735	0.0008426	0.01265	0.01287	0.007142]';
% 
% % Ratio of the prism gradient to my gradient
% Ratio_Early = Prism_Gradient_Early ./ Gradient_Base;
% 
% % differences in the early gradient
% Difference_Early = Prism_Gradient_Early - Gradient_Base;
% avg_diff_Early = mean(Difference_Early)
% 
% % paired t-test on two gradients calculated
% [h_e,p_e,ci_e,stats_e] = ttest(Prism_Gradient_Early,Gradient_Base);
% 
% 
% %% stats for late
% 
% % Slope values taken from Prism 
% Prism_Gradient_late = [0.002613	0.0009840	0.0004476	-0.006423	...
%     0.001429		-0.006423	0.001845	0.008591	-0.003985	...
%     0.02944	0.008351	-0.0009682	-0.001633	0.02360		0.02561	...
%     -0.003903	0.002709	-0.01879	0.01463	0.001429	0.01051	...
%     0.0008773	-0.001303		0.0003190	-0.001236	0.004620	...
%     -0.004705	0.01288	0.02958	0.004008	0.02511	0.007929	...
%     0.01480		-0.0003569	-0.004701	-0.004183	-0.0003709	...
%     -9.999e-005	0.01713	0.0004761	0.005086	0.002419	0.007838]';
% 
% % Ratio between the two 
% Ratio_Late = Prism_Gradient_late ./ Gradient_Post;
% 
% % Differences between the two
% Difference_late = Prism_Gradient_late - Gradient_Post;
% avg_diff_late = mean(Difference_late)
% 
% % paired t-test on the two gradients
% [h_l,p_l,ci_l,stats_l] = ttest(Prism_Gradient_late,Gradient_Post);