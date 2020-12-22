%% Gradient Plot Script
% heavily uses plotg
% (https://www.mathworks.com/matlabcentral/fileexchange/67385-plotg) and
% the support listed here
% (https://www.mathworks.com/matlabcentral/answers/265914-how-to-make-my-own-custom-colormap)

% consider making more user friendly in future if required
% define a vector of points for the colours to be present at
vec = [100; 75; 25; 0];
% define a vector of the hexidecimal colours you want at the points above
hex = ['#0072BD';'#0072BD';'#D95319';'#D95319'];
% convert to decimal
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
% interpolate to get a blended colormap
N = 128;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
% use plotg to plot your line with the map you've just defined
plotg(S.array(:,1),S.array(:,2),[],map)