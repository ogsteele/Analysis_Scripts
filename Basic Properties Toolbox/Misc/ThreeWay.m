function BarPlot3(x,c,e)

clear all

b = bar(c,a);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [0 1 0];
b.CData(3,:) = [1 0 0];
hold on
errorbar(a,error,'.');

