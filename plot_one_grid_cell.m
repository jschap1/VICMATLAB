%% 

load('FLUXES.mat')

%% Time series plot of a particular flux variable

figure
plot(FLUXES.time, FLUXES.ts.cell_37_59375_120_90625_txt.aero_resist)
xlabel('time'); ylabel('aero_resist (s/m)');
set(gca, 'FontSize', 14)  

% saveas(gcf, fullfile(saveloc, 'savename'));
% savefig(gcf, fullfile(saveloc, 'savename'));

