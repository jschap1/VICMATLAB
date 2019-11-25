% Plot water balance 
%
% Plots the water balance terms from the VIC model
% INPUTS
% info = maps and/or time series of VIC outputs
% dat = data, either maps or time series
% type = 'maps' or 'ts'

function plotwb(info, dat, type)

fontsize = 14;

if strcmp(type, 'maps')
    
    % fluxes
    precipitation = xyz2grid(info.lon, info.lat, dat.OUT_PREC);
    evaporation = xyz2grid(info.lon, info.lat, dat.OUT_EVAP);
    runoff = xyz2grid(info.lon, info.lat, dat.OUT_RUNOFF);
    baseflow = xyz2grid(info.lon, info.lat, dat.OUT_BASEFLOW);
    
    figure
    subplot(2,2,1)
    plotraster(info.lon, info.lat, precipitation, 'Precipitation (mm)', 'Lon', 'Lat')
    subplot(2,2,2)
    plotraster(info.lon, info.lat, evaporation, 'Evaporation (mm)', 'Lon', 'Lat')
    subplot(2,2,3)
    plotraster(info.lon, info.lat, runoff, 'Runoff (mm)', 'Lon', 'Lat')
    subplot(2,2,4)
    plotraster(info.lon, info.lat, baseflow, 'Baseflow (mm)', 'Lon', 'Lat')

    % states
    swe = xyz2grid(info.lon, info.lat, dat.OUT_SWE);
    sm1 = xyz2grid(info.lon, info.lat, dat.OUT_SOIL_MOIST_0);
    sm2 = xyz2grid(info.lon, info.lat, dat.OUT_SOIL_MOIST_1);
    sm3 = xyz2grid(info.lon, info.lat, dat.OUT_SOIL_MOIST_2);
    
    figure
    subplot(2,2,1)
    plotraster(info.lon, info.lat, swe, 'SWE (mm)', 'Lon', 'Lat')
    subplot(2,2,2)
    plotraster(info.lon, info.lat, sm1, 'Soil moisture (1) (mm)', 'Lon', 'Lat')
    subplot(2,2,3)
    plotraster(info.lon, info.lat, sm2, 'Soil moisture (2) (mm)', 'Lon', 'Lat')
    subplot(2,2,4)
    plotraster(info.lon, info.lat, sm3, 'Soil moisture (3) (mm)', 'Lon', 'Lat')    
    
elseif strcmp(type, 'ts')

    % fluxes
    figure
    subplot(2,2,1)
    jsplot(info.time, dat.OUT_PREC, 'Precipitation', 'Time', 'Precip. (mm)', fontsize)
    subplot(2,2,2)
    jsplot(info.time, dat.OUT_EVAP, 'Evaporation', 'Time', 'Evap. (mm)', fontsize)
    subplot(2,2,3)
    jsplot(info.time, dat.OUT_RUNOFF, 'Runoff', 'Time', 'Runoff (mm)', fontsize)
    subplot(2,2,4)
    jsplot(info.time, dat.OUT_BASEFLOW, 'Baseflow', 'Time', 'Q_{bas} (mm)', fontsize)
    
    % states
    figure
    subplot(2,2,1)
    jsplot(info.time, dat.OUT_SWE, 'SWE', 'Time', 'SWE (mm)', fontsize)
    subplot(2,2,2)
    jsplot(info.time, dat.OUT_SOIL_MOIST_0, 'Soil moisture (1)', 'Time', 'SM (mm)', fontsize)
    subplot(2,2,3)
    jsplot(info.time, dat.OUT_SOIL_MOIST_1, 'Soil moisture (2)', 'Time', 'SM (mm)', fontsize)
    subplot(2,2,4)
    jsplot(info.time, dat.OUT_SOIL_MOIST_2, 'Soil moisture (3)', 'Time', 'SM (mm)', fontsize)
    
end

return