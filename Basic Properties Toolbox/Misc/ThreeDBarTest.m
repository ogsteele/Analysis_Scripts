function ThreeDBarTest(x,Well) 

spikecounts =  [0,(Well(x).channel(1).numSpikes),(Well(x).channel(2).numSpikes),0;...
                (Well(x).channel(3).numSpikes),(Well(x).channel(4).numSpikes),(Well(x).channel(5).numSpikes),(Well(x).channel(6).numSpikes);...
                (Well(x).channel(7).numSpikes),(Well(x).channel(8).numSpikes),(Well(x).channel(9).numSpikes),(Well(x).channel(10).numSpikes);...
                0,(Well(x).channel(11).numSpikes),(Well(x).channel(12).numSpikes),0];
            
b = bar3(spikecounts);
colorbar

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

set(gcf,'color','w');

ylabel('Y Electrode Coordinate')
xlabel('X Electrode coordinate')
zlabel('Spike Frequency (Hz)')

end