function [f1, f2, f3, f4] = plot_difference_ts(OUTPUTS1, OUTPUTS2, figdir)

fontsize = 18;
height1 = 500;
width1 = 700;

% Water balance variables
% Fluxes
f1 = figure;
set(f1, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_PREC, OUTPUTS2.WB.ts.OUT_PREC, 'Precipitation', 'Time', 'Precipitation (mm)', fontsize)
grid on
subplot(2,2,2)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_EVAP, OUTPUTS2.WB.ts.OUT_EVAP, 'Evaporation', 'Time', 'Evaporation (mm)', fontsize)
grid on
subplot(2,2,3)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_RUNOFF, OUTPUTS2.WB.ts.OUT_RUNOFF, 'Runoff', 'Time', 'Runoff (mm)', fontsize)
grid on
subplot(2,2,4)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_BASEFLOW, OUTPUTS2.WB.ts.OUT_BASEFLOW, 'Baseflow', 'Time', 'Baseflow (mm)', fontsize)
grid on
saveas(f1, fullfile(figdir, 'water_balance_fluxes_diff.png'))

% Spatial average time series
% Water balance variables
% States
f2 = figure;
set(f2, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SWE, OUTPUTS2.WB.ts.OUT_SWE, 'SWE', 'Time', 'SWE (mm)', fontsize)
grid on
subplot(2,2,2)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_0, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_0, 'Soil moisture (layer 1)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
subplot(2,2,3)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_1, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_1, 'Soil moisture (layer 2)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
subplot(2,2,4)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_2, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_2, 'Soil moisture (layer 3)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
saveas(f2, fullfile(figdir, 'water_balance_states_diff.png'))

% Spatial average time series
% Energy balance variables
% Fluxes
f3 = figure;
set(f3, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_R_NET, OUTPUTS2.EB.ts.OUT_R_NET, 'Net radiation', 'Time', 'Rnet (W/m^2)', fontsize)
grid on
subplot(2,2,2)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_LATENT, OUTPUTS2.EB.ts.OUT_LATENT, 'Latent heat', 'Time', 'LE (W/m^2)', fontsize)
grid on
subplot(2,2,3)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_SENSIBLE, OUTPUTS2.EB.ts.OUT_SENSIBLE, 'Sensible heat', 'Time', 'H (W/m^2)', fontsize)
grid on
subplot(2,2,4)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_GRND_FLUX, OUTPUTS2.EB.ts.OUT_GRND_FLUX, 'Ground heat flux', 'Time', 'G (W/m^2)', fontsize)
grid on
saveas(f3, fullfile(figdir, 'energy_balance_fluxes_diff.png'))

% Spatial average time series
% Energy balance variables
% States
f4 = figure;
set(f4, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_ALBEDO, OUTPUTS2.EB.ts.OUT_ALBEDO, 'Albedo', 'Time', '\alpha (-)', fontsize)
grid on
subplot(2,2,2)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_SURF_TEMP, OUTPUTS2.EB.ts.OUT_SURF_TEMP, 'Surface temperature', 'Time', 'T_{surf} (deg. C)', fontsize)
grid on
subplot(2,2,3)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_AIR_TEMP, OUTPUTS2.EB.ts.OUT_AIR_TEMP, 'Air temperature', 'Time', 'T_{air} (deg. C)', fontsize)
grid on
subplot(2,2,4)
plotdiff_ts(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_REL_HUMID, OUTPUTS2.EB.ts.OUT_REL_HUMID, 'Relative humidity', 'Time', 'RH (%)', fontsize)
grid on
saveas(f4, fullfile(figdir, 'energy_balance_states_diff.png'))

return