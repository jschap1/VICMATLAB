% Plot energy balance 
%
% Plots the energy balance terms from the VIC model
% INPUTS
% info = maps and/or time series of VIC outputs
% dat = data, either maps or time series
% type = 'maps' or 'ts'

function ploteb(info, dat, type)

fontsize = 14;

if strcmp(type, 'maps')

    % fluxes
    Rn = xyz2grid(info.lon, info.lat, dat.OUT_R_NET);
    H = xyz2grid(info.lon, info.lat, dat.OUT_SENSIBLE);
    LE = xyz2grid(info.lon, info.lat, dat.OUT_LATENT);
    G = xyz2grid(info.lon, info.lat, dat.OUT_GRND_FLUX);
    
    figure
    subplot(2,2,1)
    plotraster(info.lon, info.lat, Rn, 'Net radiation (W/m^2)', 'Lon', 'Lat')
    subplot(2,2,2)
    plotraster(info.lon, info.lat, LE, 'Latent heat (W/m^2)', 'Lon', 'Lat')
    subplot(2,2,3)
    plotraster(info.lon, info.lat, H, 'Sensible heat (W/m^2)', 'Lon', 'Lat')
    subplot(2,2,4)
    plotraster(info.lon, info.lat, G, 'Ground heat flux (W/m^2)', 'Lon', 'Lat')

    % states
    albedo = xyz2grid(info.lon, info.lat, dat.OUT_ALBEDO);
    air_temp = xyz2grid(info.lon, info.lat, dat.OUT_AIR_TEMP);
    sfc_temp = xyz2grid(info.lon, info.lat, dat.OUT_SURF_TEMP);
    aero_resist = xyz2grid(info.lon, info.lat, dat.OUT_AERO_RESIST);
    
    figure
    subplot(2,2,1)
    plotraster(info.lon, info.lat, albedo, 'Albedo (-)', 'Lon', 'Lat')
    subplot(2,2,2)
    plotraster(info.lon, info.lat, air_temp, 'Air temperature (deg. C)', 'Lon', 'Lat')
    subplot(2,2,3)
    plotraster(info.lon, info.lat, sfc_temp, 'Surface temperature (deg. C)', 'Lon', 'Lat')
    subplot(2,2,4)
    plotraster(info.lon, info.lat, aero_resist, 'Aerodynamic resistance (s/m)', 'Lon', 'Lat')    
    
elseif strcmp(type, 'ts')

    % fluxes
    figure
    subplot(2,2,1)
    jsplot(info.time, dat.OUT_R_NET, 'Net radiation (W/m^2)', 'Time', 'Rn (W/m^2)', fontsize)
    subplot(2,2,2)
    jsplot(info.time, dat.OUT_LATENT, 'Latent heat flux (W/m^2)', 'Time', 'LE (W/m^2)', fontsize)
    subplot(2,2,3)
    jsplot(info.time, dat.OUT_SENSIBLE, 'Sensible heat flux (W/m^2)', 'Time', 'H (W/m^2)', fontsize)
    subplot(2,2,4)
    jsplot(info.time, dat.OUT_GRND_FLUX, 'Ground heat flux (W/m^2)', 'Time', 'G (W/m^2)', fontsize)
    
    % states
    figure
    subplot(2,2,1)
    jsplot(info.time, dat.OUT_ALBEDO, 'Albedo (-)', 'Time', 'Albedo (-)', fontsize)
    subplot(2,2,2)
    jsplot(info.time, dat.OUT_AIR_TEMP, 'Air temperature (deg. C)', 'Time', 'Temperature (deg. C)', fontsize)
    subplot(2,2,3)
    jsplot(info.time, dat.OUT_SURF_TEMP, 'Surface temperature (deg. C)', 'Time', 'Temperature (deg. C)', fontsize)
    subplot(2,2,4)
    jsplot(info.time, dat.OUT_AERO_RESIST, 'Aerodynamic resistance (s/m)', 'Time', 'r_a (s/m)', fontsize)
       
    
end


return