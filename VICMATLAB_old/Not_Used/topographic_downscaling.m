% Topographic downscaling of meteorological forcing data
%
% Written 5/14/2019 JRS
%
% Fits in the prepare_merra2_forcings script as follows
% for d=begin_ind:end_ind
%     % topographic downscaling
% end
%    
% variables to output
% [temperature, precipitation, pressure, shortwave, longwave, vapor_pressure, wind_speed]

resolution = target_res;

% constants
sigma_sb = 5.670e-8; % Stephen-Boltzmann constantcompu
gamma_lapse = 6.5/1000; % moist adiabatic lapse rate, K/m
Rd = 286.9; % J/kgK
Rv = 461.4; % J/kgK
g = 9.81; % gravitational acceleration m/s^2

% function handles
qv2vp = @(qv,ps) (qv.*ps)./(qv+0.622); % vp units are same as ps units
vp2qv = @(vp, ps) 0.622*vp./(ps - vp);
get_resultant = @(x,y) (x.^2 + y.^2).^(1/2);
vp2tdew = @(vp) (log(vp) + 0.49299)./(0.0707 - 0.00421.*log(vp)); % kPa/deg. C
tdew2vp = @(tdew) 0.6108*exp((17.27*tdew)./(237.3+tdew)); % kPa/deg. C

%% Temperature

% crop coarse resolution MERRA-2 temperature to study region
temp_crop = zeros(24, height, width);
return_cell = 0;
temperature = hdfread(fullfile(merra_dir, air_names(1).name), 'EOSGRID', 'Fields', 'T2M');
for h=1:24
    temp_crop(h,:,:) = imcrop(squeeze(temperature(h,:,:)), rect);
end

elev_crop = imcrop(elev_raster, rect); % elevation data from MERRA-2

% lapse temperature to a reference level (z=0)
temp_0 = zeros(size(temp_crop));
for k=1:24
    temp_0(k, :, :) = squeeze(temp_crop(k,:,:)) + gamma_lapse*elev_crop;
end

% downscale temperature to fine resolution
[temp0_fine, target_lon, target_lat] = interpolate_merra2(temp_0, target_res, out_lon, out_lat, 'linear', return_cell);
[nr, nc, ~] = size(temp0_fine);
ncells = nr*nc;

% check the downscaled T0 map
figure, 
imagesc(target_lon, target_lat, squeeze(temp0_fine(:,:,1)))
set(gca, 'ydir', 'normal')

% use fine resolution elevation data to "un-lapse" temperature
[elev_fine, R_elev_fine] = geotiffread('/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif');

R_fine = makerefmat(target_lon(1), target_lat(1), resolution, resolution);
geotiffwrite('./Data/IRB/Experimental/temperature_example.tif', squeeze(temp0_fine(:,:,1)), R_fine)

[elev_fine_crop, R_elev_fine_crop] = geotiffread('./Data/IRB/Experimental/cropped_DEM_for_downscaling_1_16.tif');
elev_fine_crop = flipud(elev_fine_crop);
elev_fine_crop(elev_fine_crop < -9999) = NaN;
figure
imagesc(target_lon, target_lat, elev_fine_crop) 
set(gca, 'ydir', 'normal')

temp_fine = zeros(size(temp0_fine));
for k=1:24
    temp_fine(:, :, k) = squeeze(temp0_fine(:,:,k)) - gamma_lapse*elev_fine_crop;
end

temp_fine = temp_fine - 273.15; % convert to deg. C

% check the topographically downscaled T map
figure, 
imagesc(target_lon, target_lat, squeeze(temp_fine(:,:,1)))
set(gca, 'ydir', 'normal'), colorbar

% Figure for prospectus

irb = shaperead('/Volumes/HD3/SWOTDA/Data/IRB/Experimental/irb_merit.shp');

%%
figure
subplot(1,3,1)
imagesc(target_lon, target_lat, squeeze(temp_crop(12,:,:)) - 273.15)
% hold on
% plot(irb(1).X, irb(1).Y, 'black', 'linewidth', 2)
% hold off
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Temperature (degrees C)'), colorbar
xlabel('Lon')
ylabel('Lat')

subplot(1,3,2)
imagesc(target_lon, target_lat, elev_fine_crop)
% hold on
% plot(irb(1).X, irb(1).Y, 'black', 'linewidth', 2)
% hold off
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Elevation (m)'), colorbar
xlabel('Lon')
ylabel('Lat')

subplot(1,3,3)
imagesc(target_lon, target_lat, temp_fine(:,:,12))
% hold on
% plot(irb(1).X, irb(1).Y, 'black', 'linewidth', 2)
% hold off
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Temperature (degrees C)'), colorbar
xlabel('Lon')
ylabel('Lat')

%% Precipitation

% Use the numbers from Liston and Elder (2006), table 1
% This method comes from Thornton et al. (1997)
% precipitation % mm/hr

% precipitation scaling factor
xi = [0.35, 0.35, 0.35, 0.3, 0.25, 0.2, 0.2, 0.2, 0.2, 0.25, 0.3, 0.35];

% crop coarse resolution MERRA-2 precipitation to study region
prec_crop = zeros(24, height, width);
for h=1:24
    prec_crop(h,:,:) = imcrop(squeeze(precipitation(h,:,:)), rect);
end

% interpolate precipitation to fine resolution
prec_interp = interpolate_merra2(prec_crop, target_res, out_lon, out_lat, 'nearest', return_cell);

% check the interpolated precipitation map
figure
imagesc(target_lon, target_lat, squeeze(prec_interp(:,:,1)))
title('interpolated precipitation')
set(gca, 'ydir', 'normal')

% check the month to get the appropriate xi factor
xi_m = xi(month(forc_dates(d)));

% interpolate MERRA-2 elevations to fine resolution
e2 = zeros(1, coarse_c, coarse_r);
e2(1,:,:) = elev_crop'; % formatting the input for interpolate_merra2
elev_interp = interpolate_merra2(e2, target_res, out_lon, out_lat, 'nearest', return_cell);

% check interpolated elevation map
figure
imagesc(target_lon, target_lat, squeeze(elev_interp))
title('interpolated elevation')
set(gca, 'ydir', 'normal')

% re-scale the precipitation according to elevation and time of year
dz = elev_fine_crop - elev_interp;
K = (1+xi_m*dz)./(1-xi_m*dz);

% this method is faster and gives equivalent results to the commented-out loop
prec_fine = K.*prec_interp; 

% prec_fine = zeros(size(prec_interp));
% for h=1:24
%     prec_fine(:,:,h) = K.*prec_interp(:,:,h);
% end

% remove negative precipitation values
prec_fine(prec_fine<0) = 0;

% apply the basinmask
mask_logical = logical(basinmask);
pp = squeeze(prec_fine(:,:,1));
pp(~mask_logical) = NaN;

% check the topographically downscaled precipitation map
figure, 
imagesc(target_lon, target_lat, pp)
set(gca, 'ydir', 'normal')

%% Vapor pressure

% vapor_pressure % kPa

% crop coarse resolution MERRA-2 vapor pressure to study region
vp_crop = zeros(24, height, width);
return_cell = 0;
for h=1:24
    vp_crop(h,:,:) = imcrop(squeeze(vapor_pressure(h,:,:)), rect);
end

% convert to dewpoint temperature (VP in kPa, Tdew in deg. C)
tdew = vp2tdew(vp_crop) + 273.15; % K

% lapse dewpoint temperature to a reference level (z=0)
tdew0 = zeros(size(tdew));
for k=1:24
    tdew0(k, :, :) = squeeze(tdew(k,:,:)) + gamma_lapse*elev_crop;
end

% downscale dewpoint temperature to fine resolution
[tdew0_fine, ~, ~] = interpolate_merra2(tdew0, target_res, out_lon, out_lat, 'linear', return_cell);

% check the downscaled T_dew map
figure, 
imagesc(target_lon, target_lat, squeeze(temp0_fine(:,:,1)))
set(gca, 'ydir', 'normal')

% use fine resolution elevation data to "un-lapse" dewpoint temperature
tdew_fine = zeros(size(tdew0_fine));
for k=1:24
    tdew_fine(:, :, k) = squeeze(tdew0_fine(:,:,k)) - gamma_lapse*elev_fine_crop;
end

% convert to vapor pressure
vp_fine = tdew2vp(tdew_fine-273.15);

% check the topographically downscaled vapor pressure map
figure, 
imagesc(target_lon, target_lat, squeeze(vp_fine(:,:,1)))
set(gca, 'ydir', 'normal')

%% Downwelling longwave (must have already downscaled vapor pressure and temperature)

eo_fine = vp_fine.*10; % surface level vapor pressure (mb), converting from kPa to mb 
emissivity_fine =  0.179*(eo_fine./100).^(1/7).*(exp(350./(temp_fine+273.15))); % temperature (K)

% check emissivity
figure, 
imagesc(target_lon, target_lat, squeeze(emissivity_fine(:,:,1)))
set(gca, 'ydir', 'normal'), colorbar
title('Atmos. emissivity')

longwave_fine = emissivity_fine.*sigma_sb.*(temp_fine+273.15).^4; % (W/m^2)

% check incoming longwave radiation
figure, 
imagesc(target_lon, target_lat, squeeze(longwave_fine(:,:,1)))
set(gca, 'ydir', 'normal'), colorbar
title('Incoming longwave (W/m^2)')

%% Pressure (must have already downscaled temperature)

% pressure % kPa
pressure = hdfread(fullfile(merra_dir, air_names(1).name), 'EOSGRID', 'Fields', 'PS');

% crop coarse resolution MERRA-2 pressure to study region
ps_crop = zeros(24, height, width);
for h=1:24
    ps_crop(h,:,:) = imcrop(squeeze(pressure(h,:,:)), rect);
end

exp1 = g/(Rd*gamma_lapse);

% lapse pressure to a reference level (z=0)
p0 = zeros(size(tdew));
for k=1:24    
    ps_hr = squeeze(ps_crop(k,:,:));
    t0_hr = squeeze(temp_0(k,:,:));

    p0(k, :, :) = ps_hr.*(t0_hr./(t0_hr - gamma_lapse*elev_crop)).^(exp1);
end

% check the reference pressure map
figure, 
imagesc(target_lon, target_lat, squeeze(p0(k, :, :)))
set(gca, 'ydir', 'normal')

% downscale pressure to fine resolution
[p0_fine, ~, ~] = interpolate_merra2(p0, target_res, out_lon, out_lat, 'linear', return_cell);

% check the downscaled P0 map
figure
imagesc(target_lon, target_lat, squeeze(p0_fine(:,:,1)))
set(gca, 'ydir', 'normal')

% use fine resolution elevation data to "un-lapse" pressure
ps_fine = zeros(size(p0_fine));
for k=1:24
    p0_hr = squeeze(p0_fine(:,:,k));
    t0_hr = squeeze(temp0_fine(:,:,k)); 
    ps_fine(:, :, k) = p0_hr.*(t0_hr./(t0_hr - gamma_lapse*elev_fine_crop)).^(-exp1);
end

% check the topographically downscaled pressure map
figure, 
imagesc(target_lon, target_lat, squeeze(ps_fine(:,:,1)))
set(gca, 'ydir', 'normal')


%% Wind

% requires slope, aspect, and curvature of the terrain
% Use './Data/IRB/Experimental/cropped_DEM_for_downscaling_1_16.tif'
% calculate the needed variables in GRASS GIS

NA_value = -9999;
slope = flipud(geotiffread('./Data/IRB/Experimental/cropped_slope_for_downscaling_1_16.tif'));
aspect = flipud(geotiffread('./Data/IRB/Experimental/cropped_aspect_for_downscaling_1_16.tif'));
curv = flipud(geotiffread('./Data/IRB/Experimental/cropped_curvature_for_downscaling_1_16.tif'));
slope(slope==NA_value) = NaN;
aspect(aspect==NA_value) = NaN;
curv(curv==NA_value) = NaN;

% Plot topography parameters
figure
subplot(2,2,1)
imagesc(target_lon, target_lat, elev_fine_crop)
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Elevation (m)'), colorbar

subplot(2,2,2)
imagesc(target_lon, target_lat, slope)
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Slope (degrees)'), colorbar

subplot(2,2,3)
imagesc(target_lon, target_lat, aspect)
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Aspect (degrees)'), colorbar

subplot(2,2,4)
imagesc(target_lon, target_lat, curv)
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Curvature (1/m)'), colorbar

% convert to radians
% slope = slope*pi/180;
% aspect = aspect*pi/180;

wind_speed_x = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'U2M');
wind_speed_y = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'V2M');

% crop coarse resolution MERRA-2 wind speeds to study region
wind_x_crop = zeros(24, height, width);
wind_y_crop = zeros(24, height, width);
for h=1:24
    wind_x_crop(h,:,:) = imcrop(squeeze(wind_speed_x(h,:,:)), rect);
    wind_y_crop(h,:,:) = imcrop(squeeze(wind_speed_y(h,:,:)), rect);
end

wind_speed = get_resultant(wind_x_crop, wind_y_crop);

% Check wind speed
figure
imagesc(target_lon, target_lat, squeeze(wind_speed(1,:,:)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Wind speed (m/s)'), colorbar

% calculate wind direction
wind_dir = atan2(wind_y_crop, wind_x_crop)*180/pi; % degrees
wind_dir(wind_dir<0) = wind_dir(wind_dir<0) + 360; % just for consistency with the slope and aspect

% Check wind direction
figure
imagesc(target_lon, target_lat, squeeze(wind_dir(1,:,:)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Wind direction (degrees)'), colorbar

% interpolate wind to fine resolution with nearest neighbor method
wind_dir_interp = interpolate_merra2(wind_dir, target_res, out_lon, out_lat, 'nearest', return_cell);
wind_speed_interp = interpolate_merra2(wind_speed, target_res, out_lon, out_lat, 'nearest', return_cell);

% Check interpolated wind speed
figure
imagesc(target_lon, target_lat, squeeze(wind_speed_interp(:,:,1)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Wind speed (m/s)'), colorbar

% calculate slope in wind direction
delta_angle = wind_dir_interp - aspect; % degrees
wind_slope = slope.*cos(delta_angle*pi/180);

% Check slope in wind direction
figure
imagesc(target_lon, target_lat, squeeze(wind_slope(:,:,1)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Slope in direction (degrees)'), colorbar
    
% scale curvature and wind slope to be between -0.5 and 0.5
scaled_curv = nma_rescale(curv,-0.5,0.5);
scaled_wind_slope = nma_rescale(wind_slope,-0.5,0.5);

% Adjust wind to account for topography
weighting_factor = 1 + 0.42*scaled_curv + 0.58*scaled_wind_slope;
wind_fine = weighting_factor.*wind_speed_interp;

% Check weighting factor
figure
imagesc(target_lon, target_lat, squeeze(weighting_factor(:,:,1)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Weighting factor'), colorbar

% Check downscaled wind speed
figure
imagesc(target_lon, target_lat, squeeze(wind_fine(:,:,1)))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
title('Wind speed (m/s)'), colorbar

%% Downwelling shortwave

shortwave = hdfread(fullfile(merra_dir, rad_names(1).name), 'EOSGRID', 'Fields', 'SWGDN');

% crop coarse resolution MERRA-2 shortwave radiation to study region
shortwave_crop = zeros(24, height, width);
for h=1:24
    shortwave_crop(h,:,:) = imcrop(squeeze(shortwave(h,:,:)), rect);
end

% do downscaling 

% simple interpolation ------------------------------------------------
[sw_fine, ~, ~] = interpolate_merra2(shortwave_crop, target_res, out_lon, out_lat, 'nearest', return_cell);

% topographic downscaling (following MOD-WET) -------------------------

% Calculate slope and aspect
% [slope, aspect] = generate_slope_and_aspect_from_DEM(elev_fine_crop, target_lon, target_lat);
% this function requires a projected coordinate system, I believe. Use
% slope and aspect calculated in GIS, instead

slope = flipud(geotiffread('./Data/IRB/Experimental/cropped_slope_for_downscaling_1_16.tif'));
slope(slope == -9999) = NaN;

aspect = flipud(geotiffread('./Data/IRB/Experimental/cropped_aspect_for_downscaling_1_16.tif'));
aspect(aspect == -9999) = NaN;
aspect = 180 - aspect - 90 + 360;
aspect(aspect > 360) = aspect(aspect > 360) - 360;

% make sure the aspect convention is the ArcGIS convention. Will need to
% convert from aspect calculated in GRASS. Slope should be in degrees.

% phi_g = 90; % GRASS GIS convention
% phi_a = 180 - phi_g - 90 + 360; % ArcMap convention

maskfile = './Data/IRB/Experimental/cropped_basinmask_for_downscaling_1_16.tif';
[mask, Rmask] = geotiffread(maskfile);
mask = flipud(mask);
mask(mask<0) = 0;

figure

subplot(1,2,1)
imagesc(target_lon, target_lat, mask.*slope)
title('Slope')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(1,2,2)
imagesc(target_lon, target_lat, mask.*aspect)
title('Aspect')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% Calculate shade and sky-view factor
[shade_lookup_table,SVF,discrete_zenith_values,discrete_azimuth_values] =  ...
    compute_shade_lookup_table_and_SVF(target_lon, target_lat, elev_fine_crop, slope, aspect);

figure
subplot(2,2,1)
imagesc(target_lon, target_lat, mask.*SVF)
title('Sky view factor')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')

subplot(2,2,2)
imagesc(target_lon, target_lat, mask.*elev_fine_crop)
title('Elevation (m)')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')

subplot(2,2,3)
imagesc(target_lon, target_lat, mask.*slope)
title('Slope (degrees)')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')

subplot(2,2,4)
imagesc(target_lon, target_lat, mask.*aspect)
title('Aspect (degrees)')
xlabel('Longitude');
ylabel('Latitude');
colorbar
set(gca, 'ydir', 'normal')

% Calculate TOA incoming flux (W/m^2)

% Calculate time zone shift for each grid cell
tz_shift = zeros(nc, 1);
tz_shift_mat = zeros(nr, nc);
for c=1:nc
    tz_shift(c) = calc_tz_shift(target_lon(c));
    tz_shift_mat(:, c) = tz_shift(c);
end

% DOY = day of year
% UTC = time in hours

DOY = day(forc_dates(d), 'dayofyear'); % make sure you don't define any variables called "day" or this won't work
% DOY = 149; % filling in some number for now

UTC = 1:24;

RsTOA = NaN(nr, nc, 24);
zenith_deg = NaN(nr, nc, 24);
azimuth_deg = NaN(nr, nc, 24);
sunrise = NaN(nr, nc);
sunset = NaN(nr, nc);
solar_decl_rad = NaN(nr, nc, 24);
hour_angle_rad = NaN(nr, nc, 24);

for r=1:nr
    for c=1:nc
        [RsTOA(r,c,:), zenith_deg(r,c,:), azimuth_deg(r,c,:), sunrise(r,c), sunset(r,c), solar_decl_rad(r,c), hour_angle_rad(r,c,:)] = ...
            TOA_incoming_solar(DOY, UTC, tz_shift_mat(r,c), target_lat(r), target_lon(c));
    end
end

% Plot RsTOA
figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, mask.*RsTOA(:,:,h))
    title(['RsTOA (hour ' num2str(h) ' UTC'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end

% Topographic shade calculations
% interpolate shade from shade lookup table

X1 = 1:nr;
X2 = 1:nc;
X3 = discrete_zenith_values';
X4 = discrete_azimuth_values';
V = double(shade_lookup_table);
Xq1 = 1:nr;
Xq2 = 1:nc;

% shade=interpn(1:size(shade_lookup_table,1),1:size(shade_lookup_table,2),...
% discrete_zenith_values',discrete_azimuth_values',...
% double(shade_lookup_table),1:size(shade_lookup_table,1),...
% 1:size(shade_lookup_table,2),zenith_deg,azimuth_deg);

% There is definitely an issue with the shade calculation. 
% If you can get this working, I think the rest of the downscaling will
% work out, too.
shade = interpn(X1, X2, X3, X4, V, Xq1, Xq2, zenith_deg(:,:,h), azimuth_deg(:,:,h));

shade(r,c,h) = interpn(1:nr, 1:nc, discrete_zenith_values', discrete_azimuth_values', V, ...
    r, c, zenith_deg(r,c,h), azimuth_deg(r,c,h));

% This loop takes on the order of 10 minutes (started at 5:57 p.m.)
shade = NaN(nr, nc, 24);
for r=1:nr
    for c=1:nc
        for h = 1:24  
            shade(r,c,h) = interpn(1:nr, 1:nc, discrete_zenith_values', discrete_azimuth_values', V, ...
                r, c, zenith_deg(r,c,h), azimuth_deg(r,c,h));
%             Xq3 = zenith_deg(r,c,h);
%             Xq4 = azimuth_deg(r,c,h);
%             shade(r,c,h) = interpn(X1, X2, X3, X4, V, Xq1, Xq2, Xq3, Xq4);
        end
    end
    disp(r)
end
% There is no way this will be viable in a loop for every day of the
% simulation (about 14000 days). It would take nearly half a year to run on one processor.
% Need to figure out how to get the same
% result without using a loop (assuming the result is as desired)

% save the shade values (they take a while to compute)
save('./Data/irb_shade_DOY_1_1980.mat', 'shade')

figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, mask.*shade(:,:,h))
    title(['Shade (hour ' num2str(h) ' UTC'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end

% OK, yeah, I don't trust the shade calculation at all. I want to do it
% myself, from scratch, to be very sure.

% Convert angles from degrees to radians
zenith_rad = zenith_deg.*pi./180;  
azimuth_rad = azimuth_deg.*pi./180;
slope_rad = slope.*pi./180;
aspect_rad = aspect.*pi./180;

% Solar altitude angle
altitude_rad = pi/2 - zenith_rad; % radians

% Disaggregate shortwave

albedo_file = './Data/IRB/Experimental/cropped_albedo_1_16.tif';
albedo = geotiffread('./Data/IRB/Experimental/albedo_exp.tif');
albedo(albedo<0) = NaN;

[sw_fine, ~, ~] = interpolate_merra2(shortwave_crop, target_res, out_lon, out_lat, 'nearest', return_cell);

addpath('/Volumes/HD3/SWOTDA/Codes')

[Rs,RsDir,RsDif] = disaggregate_SW_JS(sw_fine, ps_fine, slope_rad, aspect_rad, ...
    zenith_rad, azimuth_rad, hour_angle_rad, shade, SVF, sunrise, sunset, RsTOA, mask, albedo);

% Direct solar radiation
figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, RsDir(:,:,h))
    title(['Direct solar radiation (W/m^2) (hour ' num2str(h) ' UTC)'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end

% Diffuse solar radiation
figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, RsDif(:,:,h))
    title(['Diffuse solar radiation (W/m^2) (hour ' num2str(h) ' UTC)'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end

% Total solar radiation
figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, Rs(:,:,h))
    title(['Direct solar radiation (W/m^2) (hour ' num2str(h) ' UTC)'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end

% Total solar radiation (simple interpolation)
figure
for h=1:24
    subplot(4,6,h)
    imagesc(target_lon, target_lat, mask.*sw_fine(:,:,h))
    title(['Direct solar radiation (W/m^2) (hour ' num2str(h) ' UTC)'])
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar
    set(gca, 'ydir', 'normal')
end


[SW0,~,~] = disaggregate_SW(shortwave, Psfc0, slope_rad, aspect_rad, zenith_rad, azimuth_rad, hour_angle_rad, shade, SVF, sunrise, sunset, RsTOA,mask,albedo);

terrain_flag
delineation_flag
slope_aspect_flag
shade_calc_flag
met_file_process_flag
easting
northing
elev
outlet_coordinate
output_filename
met_input_filename
WY_to_process
dt_met
dt_interp
met_output_mat_filename
gage_elev

MOD_WET_watershed_preprocessing(terrain_flag,delineation_flag,...
    slope_aspect_flag,shade_calc_flag,met_file_process_flag,...
        easting, northing, elev, outlet_coordinate,output_filename,...
        met_input_filename,WY_to_process,dt_met,dt_interp,...
        met_output_mat_filename,gage_elev)


% check downscaled shortwave radiation
for h=1:24
    figure
    imagesc(target_lon, target_lat, squeeze(sw_fine(:,:,h)))
    set(gca, 'ydir', 'normal')
    set(gca, 'fontsize', 18)
    title(['Shortwave radiation (W/m^2)', 'hour = ', num2str(h)]), colorbar
end



% This is still not working. Just going to do simple interpolation for now.

addpath('/Volumes/HD3/SWOTDA/Codes/SW_Disagg')
addpath(genpath('/Volumes/HD3/SWOTDA/Software/MOD_WET_2017a/MOD_WET_2017a'))



[sw0, ~, ~] = interpolate_merra2(sw_crop, target_res, out_lon, out_lat, 'nearest', return_cell);



% make lat and lon rasters
lat_mat = zeros(nr, nc);
lon_mat = zeros(nr, nc);
for r=1:nr
    for c=1:nc
        lat_mat(r,c) = target_lat(r);
        lon_mat(r,c) = target_lon(c);
    end
end

DOY = day(forc_dates(d), 'dayofyear'); 

static_data = load('./Data/IRB/Experimental/static_data.mat');
discrete_azimuth_values = static_data.discrete_azimuth_values;
discrete_zenith_values = static_data.discrete_zenith_values;
shade_lookup_table = static_data.shade_lookup_table;
SVF = static_data.SVF;

sw_fine = zeros(24, nr, nc);
for h=1:24
    sw_fine(:,:,h) = disaggregate_MERRA2_SW(sw0(:,:,h), ps_fine(:,:,h), ...
        slope, aspect, mask, albedo, DOY, h, tz_shift_mat, ...
        lat_mat, lon_mat, elev_fine_crop, ...
        shade_lookup_table, SVF, discrete_zenith_values, discrete_azimuth_values);
end



