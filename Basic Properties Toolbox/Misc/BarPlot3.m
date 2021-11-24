function BarPlot3(y,ytitle,c,e)
% y = data in the form of numerical vector 
% ytitle = character string with the ylabel present
% c = categorical 
% e = error bar
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
b = bar(c,y);
hold on
box(axes1,'on');
b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1]; % Blue Bar
    b.CData(2,:) = [0 1 0]; % Green Bar
    b.CData(3,:) = [1 0 0]; % Red Bar
hold on
    errorbar(y,e,'.','Color',[0 0 0],'LineWidth',1);
    ylabel(ytitle);
set(axes1,'FontSize',28);
end

