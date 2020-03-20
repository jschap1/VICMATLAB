%% Spatial average time series
% 
% Plots spatial average time series of water balance and energy balance
% variables and saves them to figdir

function [f1,f2,f3,f4] = plot_spatial_avg_ts_multi(OUTPUTS1, OUTPUTS2, figdir)

fontsize = 18;
height1 = 500;
width1 = 700;

% Water balance variables
% Fluxes
f1 = figure;
upper1 = 10;
set(f1, 'Position',  [100, 100, 100+width1, 100+height1])

subplot(2,2,1)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_PREC, 'Precipitation', 'Time', 'Precipitation (mm)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_PREC)
grid on
ylim([0,40])

subplot(2,2,2)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_EVAP, 'Evaporation', 'Time', 'Evaporation (mm)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_EVAP)
grid on
ylim([0,upper1])

subplot(2,2,3)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_RUNOFF, 'Runoff', 'Time', 'Runoff (mm)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_RUNOFF)
grid on
ylim([0,upper1])

subplot(2,2,4)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_BASEFLOW, 'Baseflow', 'Time', 'Baseflow (mm)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_BASEFLOW)
grid on
ylim([0,upper1])

saveas(f1, fullfile(figdir, 'water_balance_fluxes.png'))

% Spatial average time series
% Water balance variables
% States
f2 = figure;
upper2 = 400;
set(f2, 'Position',  [100, 100, 100+width1, 100+height1])

subplot(2,2,1)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SWE, 'SWE', 'Time', 'SWE (mm)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_SWE)
grid on
ylim([0,upper2])

subplot(2,2,2)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_0, 'Soil moisture (layer 1)', 'Time', 'Soil moisture (mm)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_0)
grid on
ylim([0,upper2])

subplot(2,2,3)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_1, 'Soil moisture (layer 2)', 'Time', 'Soil moisture (mm)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_1)
grid on
ylim([0,upper2])

subplot(2,2,4)
jsplot(OUTPUTS1.time, OUTPUTS1.WB.ts.OUT_SOIL_MOIST_2, 'Soil moisture (layer 3)', 'Time', 'Soil moisture (mm)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.WB.ts.OUT_SOIL_MOIST_2)
grid on
ylim([0,upper2])

saveas(f2, fullfile(figdir, 'water_balance_states.png'))

% Spatial average time series
% Energy balance variables
% Fluxes
f3 = figure;
lower3 = -50;
upper3 = 250;
set(f3, 'Position',  [100, 100, 100+width1, 100+height1])

subplot(2,2,1)
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_R_NET, 'Net radiation', 'Time', 'Rnet (W/m^2)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_R_NET)
grid on
ylim([lower3, upper3])
line([OUTPUTS1.time(1),OUTPUTS1.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')

subplot(2,2,2)
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_LATENT, 'Latent heat', 'Time', 'LE (W/m^2)', fontsize)
hold on
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_LATENT)
grid on
ylim([lower3, upper3])

subplot(2,2,3)
hold on 
line([OUTPUTS1.time(1),OUTPUTS1.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_SENSIBLE, 'Sensible heat', 'Time', 'H (W/m^2)', fontsize)
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_SENSIBLE)
grid on
ylim([lower3, upper3])

subplot(2,2,4)
hold on 
line([OUTPUTS1.time(1),OUTPUTS1.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_GRND_FLUX, 'Ground heat flux', 'Time', 'G (W/m^2)', fontsize)
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_GRND_FLUX)
grid on
ylim([lower3, upper3])

hold on 
line([OUTPUTS1.time(1),OUTPUTS1.time(end)],[0,0], 'Color', 'black', 'LineStyle', '--')
saveas(f3, fullfile(figdir, 'energy_balance_fluxes.png'))

% Spatial average time series
% Energy balance variables
% States
f4 = figure;
set(f4, 'Position',  [100, 100, 100+width1, 100+height1])

subplot(2,2,1)
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_ALBEDO, 'Albedo', 'Time', '\alpha (-)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_ALBEDO)
grid on
ylim([0,1])

subplot(2,2,2)
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_SURF_TEMP, 'Surface temperature', 'Time', 'T_{surf} (deg. C)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_SURF_TEMP)
ylim([-10,30])
grid on

subplot(2,2,3)
jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_AIR_TEMP, 'Air temperature', 'Time', 'T_{air} (deg. C)', fontsize)
hold on 
plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_AIR_TEMP)
grid on
ylim([-10,30])

if exist('OUTPUTS1.EB.maps.OUT_REL_HUMID', 'var')
    subplot(2,2,4)
    jsplot(OUTPUTS1.time, OUTPUTS1.EB.ts.OUT_REL_HUMID, 'Relative humidity', 'Time', 'RH (%)', fontsize)
    hold on
    plot(OUTPUTS2.time, OUTPUTS2.EB.ts.OUT_REL_HUMID)
    grid on
    ylim([0,100])
end

saveas(f4, fullfile(figdir, 'energy_balance_states.png'))

return