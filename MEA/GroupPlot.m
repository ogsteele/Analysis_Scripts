% Data to be plotted as a bar graph
model_series = ...
    [(Differentiation.Mean_DIV.KO(1).meanSI), (Differentiation.Mean_DIV.WT(1).meanSI), ...
    (Differentiation.Mean_DIV.MIX(1).meanSI); ...
    (Differentiation.Mean_DIV.KO(2).meanSI), (Differentiation.Mean_DIV.WT(2).meanSI), ...
    (Differentiation.Mean_DIV.MIX(2).meanSI); ...
    (Differentiation.Mean_DIV.KO(3).meanSI), (Differentiation.Mean_DIV.WT(3).meanSI), ...
    (Differentiation.Mean_DIV.MIX(3).meanSI); ...
    (Differentiation.Mean_DIV.KO(4).meanSI), (Differentiation.Mean_DIV.WT(4).meanSI), ...
    (Differentiation.Mean_DIV.MIX(4).meanSI)]

%Data to be plotted as the error bars
model_error = ...
    [(Differentiation.Mean_DIV.KO(1).nansemSI), (Differentiation.Mean_DIV.WT(1).nansemSI), ...
    (Differentiation.Mean_DIV.MIX(1).nansemSI); ...
    (Differentiation.Mean_DIV.KO(2).nansemSI), (Differentiation.Mean_DIV.WT(2).nansemSI), ...
    (Differentiation.Mean_DIV.MIX(2).nansemSI); ...
    (Differentiation.Mean_DIV.KO(3).nansemSI), (Differentiation.Mean_DIV.WT(3).nansemSI), ...
    (Differentiation.Mean_DIV.MIX(3).nansemSI); ...
    (Differentiation.Mean_DIV.KO(4).nansemSI), (Differentiation.Mean_DIV.WT(4).nansemSI), ...
    (Differentiation.Mean_DIV.MIX(4).nansemSI)]
% Creating axes and the bar graph
ax = axes;
h = bar(model_series,'BarWidth',1);
% Set color for each bar face
% h(1).FaceColor = 'blue';
% h(2).FaceColor = 'yellow';
% Properties of the bar graph as required
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1 2 3 4]);
% Naming each of the bar groups
xticklabels(ax,{ 'DIV10', 'DIV12', 'DIV14', 'DIV16'});
% X and Y labels
%xlabel ('Socio Economic Status');
ylabel ('Mean Synchronicity Index');
% Creating a legend and placing it outside the bar plot
lg = legend('KO','WT','MIX','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
box off
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
