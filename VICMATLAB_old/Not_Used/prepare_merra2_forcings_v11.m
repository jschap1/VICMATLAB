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
%
% d5 - double checked units, added comments, modified conversion formulas
%
% Updated 4/15/2019 JRS to perform cropping before interpolating to
% speed up runtime.
%
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
%
% v9
% Back to using tall arrays; but this time they are actually tall, not
% wide, so they are much faster to work with
%
% Updated 4/30/2019 v10
% Now using binary inputs and outputs for speed
%
% v11
% Final working version of MERRA-2 interpolation/VIC forcing prep code
% Note: had trouble with binary forcing files for VIC, but it works when I
% convert them to ASCII format in Matlab, then load into VIC.
% Note: Do not cd into the data directory, or Matlab will become very, very slow!
%
% TODO: make it possible to not specify an extension; VIC doesn't read in
% forcing files properly if they have extensions
%
% Add topographic downscaling capability


%% Set up environment

cd('/Volumes/HD3/SWOTDA')
addpath('/Users/jschap/Desktop/MERRA2/Codes')
addpath('/Users/jschap/Documents/Codes/VICMATLAB')

delete(gcp('nocreate')) % remove any existing parallel pools
p = parpool(); % start a parallel pool

%% Inputs

lat_range = [24 37.1875];
lon_range = [66 84];

target_res = 1/16; % can change this if not enough RAM

% MERRA-2 static file with lat/lon info
static_file = '/Users/jschap/Desktop/MERRA2/Processed/static_file_MERRA2.in'; 

% Directory to save the downscaled forcing files
% forcdir = '/Users/jschap/Desktop/MERRA2/Forc_1980-2019'; 
forcdir = './Data/IRB/VIC/Forc_2009-2019'; 

% temporary location to store the downscaled, cropped forcing data as text files
% intermediate_dir = '/Users/jschap/Desktop/MERRA2/Downscaled';
intermediate_dir = '/Volumes/HD_ExFat/downscaled_cropped_forcings';
merra_dir = '/Volumes/HD_ExFat/MERRA2';

% MERRA-2 filenames
prec_names = dir(fullfile(merra_dir, 'MERRA2*flx*.hdf'));
rad_names = dir(fullfile(merra_dir, 'MERRA2*rad*.hdf'));
air_names = dir(fullfile(merra_dir, 'MERRA2*slv*.hdf'));

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

    % variables to output
%     [temperature, precipitation, pressure, shortwave, longwave, vapor_pressure, wind_speed]

% constants
sigma_sb = 5.670e-8; % Stephen-Boltzmann constant

% function handles
qv2vp = @(qv,ps) (qv.*ps)./(qv+0.622); % vp units are same as ps units
get_resultant = @(x,y) (x.^2 + y.^2).^(1/2);

% Calculate target_lon, target_lat, nr, nc
temp_crop = zeros(24, height, width);
return_cell = 0;
temperature = hdfread(fullfile(merra_dir, air_names(1).name), 'EOSGRID', 'Fields', 'T2M');
for h=1:24
    temp_crop(h,:,:) = imcrop(squeeze(temperature(h,:,:)), rect);
end
[temp_fine, target_lon, target_lat] = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
[nr, nc, ~] = size(temp_fine);
ncells = nr*nc;

% make list of the latlon values corresponding to each grid cell, and in
% the order they appear in each of the saved temperature.txt files
xyz = raster2xyz(target_lon', target_lat', ones(nr, nc));
lonlat = xyz(:,1:2);
% hopefully the order is correct; plotting will show whether it worked

% extens = '.txt';
extens = '.bin';

M = [2^8, 2^5, 2^6, 2^3, 2^3, 2^7, 2^7];
M_TEMP = M(1);
M_PREC = M(2);
M_PS = M(3);
M_SW = M(4);
M_LW = M(5);
M_VP = M(6);
M_WIND = M(7);
% M = 1000; % multiplicative factor for writing binary data
% choose a scale factor that is a power of 2
% https://en.wikipedia.org/wiki/Scale_factor_(computer_science)
% Also, take into account the range of values you're working with

%% Loop through MERRA-2 files, day by day

% for d=1:31 % 68 minutes for 1000 days on 6 workers
parfor d=1:ndays
% (Takes on the order of one day to run on my machine)

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
    
    % unit conversions
    temperature = temperature - 273.15; % Kelvin to Celsius
    pressure = pressure./1000; % Pascal to kPa
    vapor_pressure = vapor_pressure./1000; % Pascal to kPa
    precipitation = 3600.*precipitation; % kg/m2/s to mm/hr   
    
    % Fix any spurious/negative values of precipitation
    precipitation(precipitation < 1e-5) = 0;
       
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

    temp_fine = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
    prec_fine = interpolate_merra2(prec_crop, target_res, out_lon, out_lat, return_cell);
    ps_fine = interpolate_merra2(ps_crop, target_res, out_lon, out_lat, return_cell);
    shortwave_fine = interpolate_merra2(shortwave_crop, target_res, out_lon, out_lat, return_cell);
    vp_fine = interpolate_merra2(vp_crop, target_res, out_lon, out_lat, return_cell);
    wind_fine = interpolate_merra2(wind_crop, target_res, out_lon, out_lat, return_cell);        
    
    % calculate longwave radiation (Idso, 1981) 
    % assumes a cloudless atmosphere
    eo_fine = vp_fine.*10; % surface level vapor pressure (mb), converting from kPa to mb 
    emissivity_fine =  0.179*(eo_fine./100).^(1/7).*(exp(350./(temp_fine+273.15))); % temperature (K)
    longwave_fine = emissivity_fine.*sigma_sb.*(temp_fine+273.15).^4; % (W/m^2)
    
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
    
    savename_TEMP = ['temperature_', forc_date_string(d,:), extens];
    savename_PREC = ['precip_', forc_date_string(d,:), extens];
    savename_PS = ['ps_', forc_date_string(d,:), extens];
    savename_SW = ['shortwave_', forc_date_string(d,:), extens];
    savename_LW = ['longwave_', forc_date_string(d,:), extens];
    savename_VP = ['vp_', forc_date_string(d,:), extens];
    savename_WIND = ['wind_', forc_date_string(d,:), extens];
        
    % takes about _ s to write each file with fwrite
    % this could be speed up by writing just one output file per day, with
    % all seven forcing variables
    fID_TEMP = fopen(fullfile(intermediate_dir, savename_TEMP),'w');
    fID_PREC = fopen(fullfile(intermediate_dir, savename_PREC),'w');
    fID_PS = fopen(fullfile(intermediate_dir, savename_PS),'w');
    fID_SW = fopen(fullfile(intermediate_dir, savename_SW),'w');
    fID_LW = fopen(fullfile(intermediate_dir, savename_LW),'w');
    fID_VP = fopen(fullfile(intermediate_dir, savename_VP),'w');
    fID_WIND = fopen(fullfile(intermediate_dir, savename_WIND),'w');
    
    fwrite(fID_TEMP, M_TEMP*temp_fine_mat, 'short');
    fwrite(fID_PREC, M_PREC*prec_fine_mat, 'short');
    fwrite(fID_PS, M_PS*ps_fine_mat, 'short');
    fwrite(fID_SW, M_SW*shortwave_fine_mat, 'short');
    fwrite(fID_LW, M_LW*longwave_fine_mat, 'short');
    fwrite(fID_VP, M_VP*vp_fine_mat, 'short');
    fwrite(fID_WIND, M_WIND*wind_fine_mat, 'short');
    
    fclose(fID_TEMP);
    fclose(fID_PREC);
    fclose(fID_PS);
    fclose(fID_SW);
    fclose(fID_LW);
    fclose(fID_VP);
    fclose(fID_WIND);
    
end
% 5.18 s for one loop

% intermediate_dir = '/Users/jschap/Desktop/MERRA2/Downscaled';
% intermediate_dir = '/Volumes/HD_ExFat/downscaled_2';
% this produces 400 MB of data for 5 days (takes 2 min. to run)
% for 365 days, it will produce roughly 30 GB of data and take about 2.5 hr
% I can live with that

% OK, so actually, it will take 96 hours and use 1.2 TB of storage, at
% least it appears it will. Save format could probably reduce storage use.

% Binary output may be faster

%% Load the cropped, downscaled forcing data and prepare VIC input files

precision = '%3.5f';

% calculate indices ahead of time
% d1 = zeros(ndays,1);
% d2 = zeros(ndays,1);
% for d=1:ndays
%     d1(d) = 24*(d-1)+1;
%     d2(d) = 24*d;   
% end

% datetimearray(:,1) % years
% datetimearray(:,2) % months
% datetimearray(:,3) % days
years_w_data = unique(datetimearray(:,1));
nyears = length(years_w_data);
% calculate the number of days in each year (accounting for partial years)
ndays_in_year = zeros(nyears, 1);
for y=1:nyears 
    ndays_in_year(y) = sum(datetimearray(:,1) == years_w_data(y));
end


%% Reformat data and write to VIC input files

% Could split this script in two: the following code could be made 
% independent of the above sections, though it is not right now

% Choose years for which to write data
startyear = 2009;
endyear = 2012;

years2run = [find(years_w_data == startyear):find(years_w_data == endyear)];

forc_name = cell(ncells, 1);
for x=1:ncells
    forc_name{x} = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) extens]);
end

day_index = 1; % for assigning the correct filename
firstyearflag = 1;

%%
% Reading the data one year at a time
% Note: cannot use a parfor loop bc there is dependence btw successive iterations
for y=years2run 
    
    TEMP = zeros(ncells, 24*ndays_in_year(y));
    PREC = zeros(ncells, 24*ndays_in_year(y));
    PS = zeros(ncells, 24*ndays_in_year(y));
    SW = zeros(ncells, 24*ndays_in_year(y));
    LW = zeros(ncells, 24*ndays_in_year(y));
    VP = zeros(ncells, 24*ndays_in_year(y));
    WIND = zeros(ncells, 24*ndays_in_year(y));
    
    % This loop takes O(minutes) to O(hours) to get through a full year
    % Also, it uses about 36 GB of RAM
    for d=1:ndays_in_year(y)

        % this could be done outside the loop
        savename_TEMP = ['temperature_', forc_date_string(day_index,:), extens];
        savename_PREC = ['precip_', forc_date_string(day_index,:), extens];
        savename_PS = ['ps_', forc_date_string(day_index,:), extens];
        savename_SW = ['shortwave_', forc_date_string(day_index,:), extens];
        savename_LW = ['longwave_', forc_date_string(day_index,:), extens];
        savename_VP = ['vp_', forc_date_string(day_index,:), extens];
        savename_WIND = ['wind_', forc_date_string(day_index,:), extens];  
        
        day_index = day_index + 1;
        
        fID_TEMP = fopen(fullfile(intermediate_dir, savename_TEMP),'r');
        fID_PREC = fopen(fullfile(intermediate_dir, savename_PREC),'r');
        fID_PS = fopen(fullfile(intermediate_dir, savename_PS),'r');
        fID_SW = fopen(fullfile(intermediate_dir, savename_SW),'r');
        fID_LW = fopen(fullfile(intermediate_dir, savename_LW),'r');
        fID_VP = fopen(fullfile(intermediate_dir, savename_VP),'r');
        fID_WIND = fopen(fullfile(intermediate_dir, savename_WIND),'r');
        
        d1 = 24*(d-1)+1;
        d2 = 24*d;
        TEMP(:, d1:d2) = fread(fID_TEMP, [ncells, 24], 'short');
        PREC(:, d1:d2) = fread(fID_PREC, [ncells, 24], 'short');
        PS(:, d1:d2) = fread(fID_PS, [ncells, 24], 'short');
        SW(:, d1:d2) = fread(fID_SW, [ncells, 24], 'short');
        LW(:, d1:d2) = fread(fID_LW, [ncells, 24], 'short');
        VP(:, d1:d2) = fread(fID_VP, [ncells, 24], 'short');
        WIND(:, d1:d2) = fread(fID_WIND, [ncells, 24], 'short');
        
        % remove any spurious negative values
        PREC(PREC<0) = 0;
        WIND(WIND<0) = -WIND(WIND<0);
        % later, figure out how to make these values not appear in the
        % first place
        
        fclose(fID_TEMP);
        fclose(fID_PREC);
        fclose(fID_PS);
        fclose(fID_SW);
        fclose(fID_LW);
        fclose(fID_VP);
        fclose(fID_WIND);
        
    end

    % produces about 8 GB per year
    if firstyearflag
        for x=1:ncells % takes about 15 minutes
            fID_FORC = fopen(forc_name{x}, 'w');
            forcing_out = [TEMP(x,:)', PREC(x,:)', PS(x,:)', SW(x,:)', LW(x,:)', VP(x,:)', WIND(x,:)'];
            fwrite(fID_FORC, forcing_out, 'short');
            fclose(fID_FORC);            
        end
        firstyearflag = 0;
        sz_old = size(forcing_out);
    else
        % takes on the order of hours; computation time increases as the
        % files grow larger in successive years
        for x=1:ncells
            % read the previous file
            fID_FORC = fopen(forc_name{x}, 'r');
            forcing_out = fread(fID_FORC, sz_old, 'short');
            fclose(fID_FORC);
            % append the new data
            fID_FORC = fopen(forc_name{x}, 'w');
            forcing_out_appended = vertcat(forcing_out, [TEMP(x,:)', PREC(x,:)', PS(x,:)', SW(x,:)', LW(x,:)', VP(x,:)', WIND(x,:)']);
            fwrite(fID_FORC, forcing_out_appended, 'short');
            fclose(fID_FORC);
        end
        sz_old = size(forcing_out_appended);
    end

end


%% Convert binary files to ASCII files

ascii_forc_dir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/Forc_2009-2012_ascii';

% parpool()

parfor x=1:ncells
    
    fID_FORC = fopen(forc_name{x}, 'r');
    forcing_in = fread(fID_FORC, sz_old, 'short');
    fclose(fID_FORC);
    
    temp_in = forcing_in(:,1)/M_TEMP;
    prec_in = forcing_in(:,2)/M_PREC;
    ps_in = forcing_in(:,3)/M_PS;
    sw_in = forcing_in(:,4)/M_SW;
    lw_in = forcing_in(:,5)/M_LW;
    vp_in = forcing_in(:,6)/M_VP;
    wind_in = forcing_in(:,7)/M_WIND;
    
    forcing_out = [temp_in, prec_in, ps_in, sw_in, lw_in, vp_in, wind_in];
        
    ascii_forc_savename = fullfile(ascii_forc_dir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision)]);
    fID = fopen(ascii_forc_savename, 'w');
    fprintf(fID, '%f %f %f %f %f %f %f \n', forcing_out');
    fclose(fID);

end


%% Plot forcing data for the last grid cell (as a check)

% x=10000;
fID_FORC = fopen(forc_name{x}, 'r');
forcing_out = fread(fID_FORC, sz_old, 'short');
fclose(fID_FORC);

% multiply by scale factor

figure

subplot(2,4,1)
plot(forcing_out(:,1)/M_TEMP), title('Temperature'), xlabel('Time'), ylabel('Degrees C')

subplot(2,4,2)
plot(forcing_out(:,2)/M_PREC), title('Precipitation'), xlabel('Time'), ylabel('mm/hr')

subplot(2,4,3)
plot(forcing_out(:,3)/M_PS), title('Pressure'), xlabel('Time'), ylabel('kPa')

subplot(2,4,4)
plot(forcing_out(:,4)/M_SW), title('Shortwave'), xlabel('Time'), ylabel('W/m^2')

subplot(2,4,5)
plot(forcing_out(:,5)/M_LW), title('Longwave'), xlabel('Time'), ylabel('W/m^2')

subplot(2,4,6)
plot(forcing_out(:,6)/M_VP), title('Vapor Pressure'), xlabel('Time'), ylabel('kPa')

subplot(2,4,7)
plot(forcing_out(:,7)/M_WIND), title('Wind Speed'), xlabel('Time'), ylabel('m/s')

%% END OF MERRA-2 DOWNSCALING CODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
