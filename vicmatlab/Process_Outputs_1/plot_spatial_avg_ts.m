%% Spatial average time series
% 
% Plots spatial average time series of water balance and energy balance
% variables and saves them to figdir

function [f1,f2,f3,f4] = plot_spatial_avg_ts(OUTPUTS, figdir)

fontsize = 18;
height1 = 500;
width1 = 700;

% Water balance variables
% Fluxes
f1 = figure;
upper1 = 40;
set(f1, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_PREC, 'Precipitation', 'Time', 'Precipitation (mm)', fontsize)
grid on
ylim([0,125])
subplot(2,2,2)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_EVAP, 'Evaporation', 'Time', 'Evaporation (mm)', fontsize)
grid on
ylim([0,upper1])
subplot(2,2,3)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_RUNOFF, 'Runoff', 'Time', 'Runoff (mm)', fontsize)
grid on
ylim([0,upper1])
subplot(2,2,4)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_BASEFLOW, 'Baseflow', 'Time', 'Baseflow (mm)', fontsize)
grid on
ylim([0,upper1])
saveas(f1, fullfile(figdir, 'water_balance_fluxes.png'))

% Spatial average time series
% Water balance variables
% States
f2 = figure;
upper2 = 200;
set(f2, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_SWE, 'SWE', 'Time', 'SWE (mm)', fontsize)
grid on
ylim([0,3000])
subplot(2,2,2)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_SOIL_MOIST_0, 'Soil moisture (layer 1)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
ylim([0,upper2])
subplot(2,2,3)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_SOIL_MOIST_1, 'Soil moisture (layer 2)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
ylim([0,upper2])
subplot(2,2,4)
jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_SOIL_MOIST_2, 'Soil moisture (layer 3)', 'Time', 'Soil moisture (mm)', fontsize)
grid on
ylim([0,850])
saveas(f2, fullfile(figdir, 'water_balance_states.png'))

% Spatial average time series
% Energy balance variables
% Fluxes
f3 = figure;
lower3 = -50;
upper3 = 250;
set(f3, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_R_NET, 'Net radiation', 'Time', 'Rnet (W/m^2)', fontsize)
grid on
ylim([lower3, upper3])
hold on 
line([OUTPUTS.time(1),OUTPUTS.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
subplot(2,2,2)
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_LATENT, 'Latent heat', 'Time', 'LE (W/m^2)', fontsize)
grid on
ylim([lower3, upper3])
subplot(2,2,3)
hold on 
line([OUTPUTS.time(1),OUTPUTS.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_SENSIBLE, 'Sensible heat', 'Time', 'H (W/m^2)', fontsize)
grid on
ylim([lower3, upper3])
subplot(2,2,4)
hold on 
line([OUTPUTS.time(1),OUTPUTS.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_GRND_FLUX, 'Ground heat flux', 'Time', 'G (W/m^2)', fontsize)
grid on
ylim([lower3, upper3])
hold on 
line([OUTPUTS.time(1),OUTPUTS.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
saveas(f3, fullfile(figdir, 'energy_balance_fluxes.png'))

% Spatial average time series
% Energy balance variables
% States
f4 = figure;
set(f4, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_ALBEDO, 'Albedo', 'Time', '\alpha (-)', fontsize)
grid on
ylim([0,1])
subplot(2,2,2)
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_SURF_TEMP, 'Surface temperature', 'Time', 'T_{surf} (deg. C)', fontsize)
ylim([-20,30])
grid on
subplot(2,2,3)
jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_AIR_TEMP, 'Air temperature', 'Time', 'T_{air} (deg. C)', fontsize)
grid on
ylim([-20,30])

if exist('OUTPUTS.EB.maps.OUT_REL_HUMID', 'var')
    subplot(2,2,4)
    jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_REL_HUMID, 'Relative humidity', 'Time', 'RH (%)', fontsize)
    grid on
    ylim([0,100])
end

saveas(f4, fullfile(figdir, 'energy_balance_states.png'))

return