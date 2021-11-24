function BarPlots(varargin)

a = length(varargin);
x = linspace(1,a,a);

for i = 1:a
    y(i) = varargin{i};
end
    
figure1 = figure;
bar(x,y)


end

% use varargin for a variable number of input arguements