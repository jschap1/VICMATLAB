%% Spatial area-average maps
% 
% Plots spatial average time series of water balance and energy balance
% variables and saves them to figdir

function [f5,f6,f7,f8] = plot_time_avg_maps(OUTPUTS, figdir)

fontsize = 18;
height1 = 500;
width1 = 700;

% Time average maps
% Water balance variables
% Fluxes
f5 = figure;
set(f5, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_PREC, 'Precipitation (mm)', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_EVAP, 'Evaporation (mm)', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_RUNOFF, 'Runoff (mm)', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_BASEFLOW, 'Baseflow (mm)', 'Lon', 'Lat')
saveas(f5, fullfile(figdir, 'water_balance_fluxes_maps.png'))

% Time average maps
% Water balance variables
% States
f6 = figure;
set(f6, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_SWE, 'SWE (mm)', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_SOIL_MOIST_0, 'Soil moisture (layer 1) (mm)', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_SOIL_MOIST_1, 'Soil moisture (layer 2) (mm)', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_SOIL_MOIST_2, 'Soil moisture (layer 3) (mm)', 'Lon', 'Lat')
saveas(f6, fullfile(figdir, 'water_balance_states_maps.png'))

% Time average maps
% Energy balance variables
% Fluxes
f7 = figure;
set(f7, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_R_NET, 'Net radiation (W/m^2)', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_LATENT, 'Latent heat (W/m^2)', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_SENSIBLE, 'Sensible heat (W/m^2)', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_GRND_FLUX, 'Ground heat flux (W/m^2)', 'Lon', 'Lat')
saveas(f7, fullfile(figdir, 'energy_balance_fluxes_maps.png'))

% Time average maps
% Energy balance variables
% States
f8 = figure;
set(f8, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_ALBEDO, 'Albedo (-)', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_SURF_TEMP, 'Surface temperature (deg. C)', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_AIR_TEMP, 'Air temperature (deg. C)', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_REL_HUMID, 'Relative humidity (%)', 'Lon', 'Lat')
saveas(f8, fullfile(figdir, 'energy_balance_states_maps.png'))

return