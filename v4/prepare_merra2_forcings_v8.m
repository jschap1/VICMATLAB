% Prepare MERRA-2 Forcing Data
%
% Prepares forcing data for the VIC model using NetCDF MERRA-2 data
% Downscales MERRA-2 data and crops it to the study domain
%
% INPUTS
% MERRA-2 atmospheric reanalysis data
% Specifically, flx, rad, and slv files
% Fluxes: M2T1NXFLX: MERRA-2 tavg1_2d_flx_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics V5.12.4
% Radiation: M2T1NXRAD: MERRA-2 tavg1_2d_rad_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics V5.12.4
% Single Level Diagnostics: M2T1NXSLV: MERRA-2 tavg1_2d_slv_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4
%
% Updated 3/16/2019 JRS, GK
% d3 - added cropping functionality so it only writes out data for the basin extent
% d5 - double checked units, added comments, modified conversion formulas
% Updated 4/15/2019 JRS to perform cropping before interpolating to
% speed up runtime.
% Updated 4/24/2019 v4 JRS to dramatically increase speed by changing the way
% it writes out the forcing input files
%
% Updated 4/24/2019 v5
% To reduce RAM-pressure, it runs over one variable at a time
%
% Updated 4/25/2019 v6 
% Uses tall arrays to avoid having to load in all the data at once
%
% Updated 4/29/2019 v7
% Using plain old datastores, instead of tall arrays. This seems the best
% solution.
%
% v8
% Added back in all the forcing variables, was previously just working with
% temperature, for simplicity
% Removed unused code

%% Set working directory

cd('/Users/jschap/Desktop/MERRA2/Raw/Sample2')
addpath('/Users/jschap/Desktop/MERRA2/Codes')
addpath('/Users/jschap/Documents/Codes/VICMATLAB')

%% Inputs

delete(gcp('nocreate')) % remove any existing parallel pools
% p = parpool(); % start a parallel pool

lat_range = [24 37.125];
lon_range = [66 82.875];

target_res = 1/16; % can change this if not enough RAM

% MERRA-2 static file with lat/lon info
static_file = '/Users/jschap/Desktop/MERRA2/Processed/static_file_MERRA2.in'; 

% Directory to save the downscaled forcing files
forcdir = '/Users/jschap/Desktop/MERRA2/Forc_1980-2019'; 

% MERRA-2 filenames
prec_names = dir('MERRA2*flx*.hdf');
rad_names = dir('MERRA2*rad*.hdf');
air_names = dir('MERRA2*slv*.hdf');

ndays = length(prec_names);

xres = 5/8;
yres = 1/2;

%% load static data
A = load(static_file);
lonlat = [A(:,3), A(:,2)];
lon = sort(unique(lonlat(:,1)));
lat = sort(unique(lonlat(:,2)));

%% Make cropping rectangle

% must use grid cell centers for makerefmat
R1 = makerefmat(min(lon), min(lat), xres, yres);

%     this formulation guarantees whole number indices
%     [ymin, xmin] = latlon2pix(R1, min(lat)+yres*10, min(lon)+xres*5) % lat, lon

% find appropriate minimum values for the cropping rectangle
minval_opt_1 = 10;
minval_opt_2 = 10;
for p=1:400
    tmp1 = abs(min(lon)+ p*xres - lon_range(1));
    tmp2 = abs(min(lat)+ p*yres - lat_range(1));
    if tmp1 < minval_opt_1
        minval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < minval_opt_2
        minval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
minlon = min(lon)+ p_opt_1*xres;
minlat = min(lat)+ p_opt_2*yres;

% make a 1-pixel border just to be safe
minlon = minlon - xres;
minlat = minlat - yres;

% find appropriate maximum values for the cropping rectangle
maxval_opt_1 = 10;
maxval_opt_2 = 10;
for p=1:400
    tmp1 = abs(minlon + p*xres - lon_range(2));
    tmp2 = abs(minlat + p*yres - lat_range(2));
    if tmp1 < maxval_opt_1
        maxval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < maxval_opt_2
        maxval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
maxlon = minlon + p_opt_1*xres;
maxlat = minlat + p_opt_2*yres;  

% make a 1-pixel border just to be safe
maxlon = maxlon + xres;
maxlat = maxlat + yres;

[ymin, xmin] = latlon2pix(R1, minlat, minlon);
[ymax, xmax] = latlon2pix(R1, maxlat, maxlon);
width = xmax - xmin;
height = ymax - ymin;

%     rect = [floor(xmin), floor(ymin), ceil(width), ceil(height)];
rect = [xmin, ymin, width, height];

% I have to do this to get the code to work...
height = height+1;
width = width+1;

out_lat = minlat:yres:maxlat;
out_lon = minlon:xres:maxlon;
    
%% get datetime information for the MERRA-2 data
datetimearray = zeros(ndays, 4);
for d=1:ndays
    air_name = air_names(d).name;
    mn_split = strsplit(air_name, '.');
    mn_date = mn_split{3};
    yr = str2double(mn_date(1:4));
    mon = str2double(mn_date(5:6));
    day = str2double(mn_date(7:8));
    datetimearray(d,:) = [yr, mon, day, 0];
end
forc_dates = datetime(datetimearray(:,1), datetimearray(:,2), datetimearray(:,3));
forc_date_string = datestr(forc_dates); % useful for making filenames

%% Crop, downscale the forcing data (temperature)

% prec_fine = cell(ndays, 1);
% ps_fine = cell(ndays, 1);
% swdn_fine = cell(ndays, 1);
% lwdn_fine = cell(ndays, 1);
% vp_fine = cell(ndays, 1);
% wind_fine = cell(ndays, 1);

% temporary location to store the downscaled, cropped forcing data as text files
intermediate_dir = '/Users/jschap/Desktop/MERRA2/Downscaled';
merra_dir = '/Users/jschap/Desktop/MERRA2/Raw/Sample2';

% constants
sigma_sb = 5.670e-8; % Stephen-Boltzmann constant

% function handles
qv2vp = @(qv,ps) (qv.*ps)./(qv+0.622); % vp units are same as ps units
get_resultant = @(x,y) (x.^2 + y.^2).^(1/2);

return_cell = 0;
for d=1:ndays
    
    temperature = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'T2M');
    precipitation = hdfread(fullfile(merra_dir, prec_names(d).name), 'EOSGRID', 'Fields', 'PRECTOT');
    pressure = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'PS');
    shortwave = hdfread(fullfile(merra_dir, rad_names(d).name), 'EOSGRID', 'Fields', 'SWGDN');
    specific_humidity = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'QV2M');
    wind_speed_x = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'U2M');
    wind_speed_y = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'V2M');

    % calculate vapor pressure
    vapor_pressure = qv2vp(specific_humidity, pressure);
    
    % calculate wind speed
    wind_speed = get_resultant(wind_speed_x, wind_speed_y); % get resultant wind vector
    
    % variables to output
%     [temperature, precipitation, pressure, shortwave, longwave, vapor_pressure, wind_speed]
    
    temp_crop = zeros(24, height, width);
    prec_crop = zeros(24, height, width);
    ps_crop = zeros(24, height, width);
    shortwave_crop = zeros(24, height, width);
    vp_crop = zeros(24, height, width);
    wind_crop = zeros(24, height, width);
    for h=1:24
        temp_crop(h,:,:) = imcrop(squeeze(temperature(h,:,:)), rect);
        prec_crop(h,:,:) = imcrop(squeeze(precipitation(h,:,:)), rect);
        ps_crop(h,:,:) = imcrop(squeeze(pressure(h,:,:)), rect);
        shortwave_crop(h,:,:) = imcrop(squeeze(shortwave(h,:,:)), rect);
        vp_crop(h,:,:) = imcrop(squeeze(vapor_pressure(h,:,:)), rect);
        wind_crop(h,:,:) = imcrop(squeeze(wind_speed(h,:,:)), rect);
    end

    if d==1
        [temp_fine, target_lon, target_lat] = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
        prec_fine = interpolate_merra2(prec_crop, target_res, out_lon, out_lat, return_cell);
        ps_fine = interpolate_merra2(ps_crop, target_res, out_lon, out_lat, return_cell);
        shortwave_fine = interpolate_merra2(shortwave_crop, target_res, out_lon, out_lat, return_cell);
        vp_fine = interpolate_merra2(vp_crop, target_res, out_lon, out_lat, return_cell);
        wind_fine = interpolate_merra2(wind_crop, target_res, out_lon, out_lat, return_cell);
        [nr, nc, ~] = size(temp_fine);
    else
        temp_fine = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
        prec_fine = interpolate_merra2(prec_crop, target_res, out_lon, out_lat, return_cell);
        ps_fine = interpolate_merra2(ps_crop, target_res, out_lon, out_lat, return_cell);
        shortwave_fine = interpolate_merra2(shortwave_crop, target_res, out_lon, out_lat, return_cell);
        vp_fine = interpolate_merra2(vp_crop, target_res, out_lon, out_lat, return_cell);
        wind_fine = interpolate_merra2(wind_crop, target_res, out_lon, out_lat, return_cell);        
    end
    
    % calculate longwave radiation (Idso, 1981)
    % assumes a cloudless atmosphere
    eo_fine = vp_fine./100; % surface level vapor pressure (mb), converting from Pascals to millibars 
    emissivity_fine =  0.179*(eo_fine./100).^(1/7).*(exp(350./temp_fine)); % temperature (K)
    longwave_fine = emissivity_fine.*sigma_sb.*temp_fine.^4; % (W/m^2)
    
    % Here, we write out temp_fine arrays for each day.
    % We will then read them back in as a
    % datastore in order to manipulate them and write VIC input files
    
    % convert to 2D
    temp_fine_mat = reshape(temp_fine, nr*nc, 24);
    prec_fine_mat = reshape(prec_fine, nr*nc, 24);
    ps_fine_mat = reshape(ps_fine, nr*nc, 24);
    shortwave_fine_mat = reshape(shortwave_fine, nr*nc, 24);
    longwave_fine_mat = reshape(longwave_fine, nr*nc, 24);
    vp_fine_mat = reshape(vp_fine, nr*nc, 24);
    wind_fine_mat = reshape(wind_fine, nr*nc, 24);
    
    savename{1} = ['temperature_', forc_date_string(d,:), '.txt'];
    savename{2} = ['precip_', forc_date_string(d,:), '.txt'];
    savename{3} = ['ps_', forc_date_string(d,:), '.txt'];
    savename{4} = ['shortwave_', forc_date_string(d,:), '.txt'];
    savename{5} = ['longwave_', forc_date_string(d,:), '.txt'];
    savename{6} = ['vp_', forc_date_string(d,:), '.txt'];
    savename{7} = ['wind_', forc_date_string(d,:), '.txt'];
    
    dlmwrite(fullfile(intermediate_dir, savename{1}), temp_fine_mat) % ncells by nt
    dlmwrite(fullfile(intermediate_dir, savename{2}), prec_fine_mat)
    dlmwrite(fullfile(intermediate_dir, savename{3}), ps_fine_mat)
    dlmwrite(fullfile(intermediate_dir, savename{4}), shortwave_fine_mat)
    dlmwrite(fullfile(intermediate_dir, savename{5}), longwave_fine_mat)
    dlmwrite(fullfile(intermediate_dir, savename{6}), vp_fine_mat)
    dlmwrite(fullfile(intermediate_dir, savename{7}), wind_fine_mat)
    
end

% this produces 400 MB of data for 5 days (takes 2 min. to run)
% for 365 days, it will produce roughly 30 GB of data and take about 2.5 hr
% I can live with that

%% Load the cropped, downscaled forcing data and prepare VIC input files

ds_temp = tabularTextDatastore(fullfile(intermediate_dir, 'temperature*'), 'FileExtension', '.txt');
ds_prec = tabularTextDatastore(fullfile(intermediate_dir, 'precip_*'), 'FileExtension', '.txt');
ds_ps = tabularTextDatastore(fullfile(intermediate_dir, 'ps_*'), 'FileExtension', '.txt');
ds_sw = tabularTextDatastore(fullfile(intermediate_dir, 'shortwave_*'), 'FileExtension', '.txt');
ds_lw = tabularTextDatastore(fullfile(intermediate_dir, 'longwave_*'), 'FileExtension', '.txt');
ds_vp = tabularTextDatastore(fullfile(intermediate_dir, 'vp_*'), 'FileExtension', '.txt');
ds_wind = tabularTextDatastore(fullfile(intermediate_dir, 'wind_*'), 'FileExtension', '.txt');

% make list of the latlon values corresponding to each grid cell, and in
% the order they appear in each of the saved temperature.txt files
xyz = raster2xyz(target_lon', target_lat', ones(nr, nc));
lonlat = xyz(:,1:2);
% hopefully the order is correct; plotting will show whether it worked

ncells = nr*nc; % number of cells in the cropped, downscaled forcing data
precision = '%3.5f';

% Do for all grid cells

% calculate indices ahead of time
d1 = zeros(ndays,1);
d2 = zeros(ndays,1);
for d=1:ndays
    d1(d) = 24*(d-1)+1;
    d2(d) = 24*d;   
end
        
% This will take a couple of days for the full domain
% A parfor loop might be helpful for cutting execution time 
tic; 
ds_temp.ReadSize = 'file';
ds_prec.ReadSize = 'file';
ds_ps.ReadSize = 'file';
ds_sw.ReadSize = 'file';
ds_lw.ReadSize = 'file';
ds_vp.ReadSize = 'file';
ds_wind.ReadSize = 'file';

%%

% trying this without datastores
% temp_names = dir(fullfile(intermediate_dir, 'temperature*'));
% prec_names = dir(fullfile(intermediate_dir, 'precip*'));
% ps_names = dir(fullfile(intermediate_dir, 'ps*'));
% shortwave_names = dir(fullfile(intermediate_dir, 'shortwave*'));
% longwave_names = dir(fullfile(intermediate_dir, 'longwave*'));
% vp_names = dir(fullfile(intermediate_dir, 'vp*'));
% wind_names = dir(fullfile(intermediate_dir, 'wind*'));

ta_temp = tall(ds_temp);

ncells = nr*nc;

idx = 1:ncells:ndays*ncells;
% t_sub = ta_temp(idx, :);
% t_sub_gather = gather(t_sub);
% t_sub_reshape = reshape(table2array(t_sub_gather)', 1, 24*ndays);
% plot(t_sub_reshape)

tic
for x=1:10
    savename1 = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) '.txt']);
    forcing_out = zeros(ndays*24,7);
    t_sub = ta_temp(idx+(x-1), :);
    t_sub_gather = gather(t_sub);
    t_sub_reshape = reshape(table2array(t_sub_gather)', 24*ndays, 1);    
    forcing_out(:,1) = t_sub_reshape;
    dlmwrite(savename1, forcing_out)
end
toc
% 3.35 seconds per cell
% seconds per cell to hours per IRB: multiply by (72960/3600)*7/6 

ncells = 10;
% for x=1:ncells
for x=1:ncells
    savename1 = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) '.txt']);
    forcing_out = zeros(ndays*24,7);
    
    reset(ds_temp);
    reset(ds_prec);
    reset(ds_ps);
    reset(ds_sw);
    reset(ds_lw);
    reset(ds_vp);
    reset(ds_wind);
    
    % read everything in before the loop, eliminate the d loop
    
    for d=1:ndays
        TEMP = read(ds_temp);
        PREC = read(ds_prec);
        PS = read(ds_ps);
        SW = read(ds_sw);
        LW = read(ds_lw);
        VP = read(ds_vp);
        WIND = read(ds_wind);
        
        % unit conversions
        TEMP = table2array(TEMP) - 273.15; % Kelvin to Celsius
        PS = table2array(PS)./1000; % Pascal to kPa
        VP = table2array(VP)./1000; % Pascal to kPa
        PREC = 3600.*table2array(PREC); % kg/m2/s to mm/hr        
        
        forcing_out(d1(d):d2(d),1) = TEMP(x,:)'; % temperature
        forcing_out(d1(d):d2(d),2) = PREC(x,:)'; % precip
        forcing_out(d1(d):d2(d),3) = PS(x,:)'; % air pressure
        forcing_out(d1(d):d2(d),4) = table2array(SW(x,:))'; % shortwave
        forcing_out(d1(d):d2(d),5) = table2array(LW(x,:))'; % longwave
        forcing_out(d1(d):d2(d),6) = VP(x,:)'; % vapor pressure
        forcing_out(d1(d):d2(d),7) = table2array(WIND(x,:))'; % wind speed
    end
    dlmwrite(savename1, forcing_out)
end
toc;

% takes about 16 seconds per grid cell per processor
% 474.7 seconds to do 100 grid cells using six workers
% 4.74 seconds per grid cell
% For 72960 grid cells, it will take 96 hours (4 days) but this is only for
% d=5. It should increase with the square of d...
% This may not be computationally feasible, either. The nested loops make
% it slow.