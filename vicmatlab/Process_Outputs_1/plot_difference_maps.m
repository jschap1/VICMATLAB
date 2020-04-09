% Plot difference maps
%
% percdiff is a flag for calculating percent differences, rather than
% absolute

function [f1, f2, f3, f4] = plot_difference_maps(OUTPUTS1, OUTPUTS2, figdir, percdiff)

height1 = 500;
width1 = 700;

% Time average maps
% Water balance variables
% Fluxes
f5 = figure;
set(f5, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_PREC, OUTPUTS2.WB.maps.OUT_PREC, 'Precipitation (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,2)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_EVAP, OUTPUTS2.WB.maps.OUT_EVAP, 'Evaporation (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,3)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_RUNOFF, OUTPUTS2.WB.maps.OUT_RUNOFF, 'Runoff (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,4)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_BASEFLOW, OUTPUTS2.WB.maps.OUT_BASEFLOW, 'Baseflow (mm)', 'Lon', 'Lat', percdiff)
saveas(f5, fullfile(figdir, 'water_balance_fluxes_diffmaps.png'))

% Time average maps
% Water balance variables
% States
f6 = figure;
set(f6, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_SWE, OUTPUTS2.WB.maps.OUT_SWE, 'SWE (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,2)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_SOIL_MOIST_0, OUTPUTS2.WB.maps.OUT_SOIL_MOIST_0, 'Soil moisture (layer 1) (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,3)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_SOIL_MOIST_1, OUTPUTS2.WB.maps.OUT_SOIL_MOIST_1, 'Soil moisture (layer 2) (mm)', 'Lon', 'Lat', percdiff)
subplot(2,2,4)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.WB.maps.OUT_SOIL_MOIST_2, OUTPUTS2.WB.maps.OUT_SOIL_MOIST_2, 'Soil moisture (layer 3) (mm)', 'Lon', 'Lat', percdiff)
saveas(f6, fullfile(figdir, 'water_balance_states_diffmaps.png'))

% Time average maps
% Energy balance variables
% Fluxes
f7 = figure;
set(f7, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_R_NET, OUTPUTS2.EB.maps.OUT_R_NET, 'Net radiation (W/m^2)', 'Lon', 'Lat', percdiff)
subplot(2,2,2)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_LATENT, OUTPUTS2.EB.maps.OUT_LATENT, 'Latent heat (W/m^2)', 'Lon', 'Lat', percdiff)
subplot(2,2,3)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_SENSIBLE, OUTPUTS2.EB.maps.OUT_SENSIBLE, 'Sensible heat (W/m^2)', 'Lon', 'Lat', percdiff)
subplot(2,2,4)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_GRND_FLUX, OUTPUTS2.EB.maps.OUT_GRND_FLUX, 'Ground heat flux (W/m^2)', 'Lon', 'Lat', percdiff)
saveas(f7, fullfile(figdir, 'energy_balance_fluxes_diffmaps.png'))

% Time average maps
% Energy balance variables
% States
f8 = figure;
set(f8, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_ALBEDO, OUTPUTS2.EB.maps.OUT_ALBEDO, 'Albedo (-)', 'Lon', 'Lat', percdiff)
subplot(2,2,2)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_SURF_TEMP, OUTPUTS2.EB.maps.OUT_SURF_TEMP, 'Surface temperature (deg. C)', 'Lon', 'Lat', percdiff)
subplot(2,2,3)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_AIR_TEMP, OUTPUTS2.EB.maps.OUT_AIR_TEMP, 'Air temperature (deg. C)', 'Lon', 'Lat', percdiff)
subplot(2,2,4)
plotdiff_map(OUTPUTS1.lon, OUTPUTS1.lat, OUTPUTS1.EB.maps.OUT_REL_HUMID, OUTPUTS2.EB.maps.OUT_REL_HUMID, 'Relative humidity (%)', 'Lon', 'Lat', percdiff)
saveas(f8, fullfile(figdir, 'energy_balance_states_diffmaps.png'))

return