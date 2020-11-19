load 'longrun/routresults.mat'

Y = routresults(:,1);
M = routresults(:,2);
D = ones(size(Y)); % data are monthly, so D = 1 is arbitrary
x = datetime(Y,M,D);

titletext = 'Monthly discharge at the basin outlet (cfs)';
saveloc = '/Users/jschapMac/Desktop/vis_Stehekin/longrun';
saveflag = 1;

figure, plot(x,routresults(:,3))

xlabel('Time')
ylabel('Flow (cfs)')
title(titletext)

if saveflag
    saveas(gcf, fullfile(saveloc, 'OUTLET_HYDROGRAPH.png'));
    savefig(gcf, fullfile(saveloc, 'OUTLET_HYDROGRAPH.fig'));
end
