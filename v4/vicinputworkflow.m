function [] = vicinputworkflow()

% INPUTS
%
% Basin mask, as list of coordinates
%
% Met. forcings at daily resolution, saved as separate NetCDF files for
% each variable and year (e.g. prec.1965.nc, wind.1990.nc)
% Basin mask with using appropriate resolution, geographic coordinates
% Soil parameter file covering >= the domain of interest
% Designed for use with Livneh met. forcing data. 
% Updated 2/12/2019 JRS
% 
% OUTPUTS
% Met. forcing file, clipped to basin boundaries
% Soil parameter file, clipped to basin boundaries
%
% Generic workflow for creating VIC input files
%
% Clips daily forcing and soil parameter data for CONUS to just the cells within 
% a specified basin shapefile. Saves the subsetted forcing and soil
% parameter data in an appropriate format to use for VIC input.
%
% Is set up to use Livneh forcing data, which is on a 444 by 922 lat/lon
% grid (1/16 degree resolution)
%
% Run the R script routinputworkflow first to generate the basin mask
% After running this script, you can run VIC4 as a disaggregator to get subdaily
% forcing inputs for VIC5.

%% Specify inputs

% Directory of daily CONUS met. forcing file
forcingpath = '/Volumes/HD3/Livneh_2013/MetNC';

% Name of basin mask at VIC modeling resolution
% basinmask = './Data/UMRB/Basin/basin.wgs.asc';
% basinmask = './Data/UMRB/rout/umrb.fract';

% Directory where clipped forcing files should be saved
forcingsavedir = '/Users/jschap/Documents/Research/Glaciers/Skagit/forc_1979-2011';

% Number of forcings in the daily CONUS met. forcing file
numforcings = 4;

% Beginning and ending years of simulation (must be included in the daily CONUS
% met. forcing file)
beginyear = 1979;
endyear = 2011;

% Number of decimal points of precision to use for forcing file names
grid_decimal = 5;

% Directory of CONUS soil parameter file
soilpath = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/';
soilname = 'soils_3L_MERIT.txt'; 

% Directory where clipped soil parameter file should be saved
soilsavedir = '.';

%% Load the mask

addpath(forcingpath)
metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');

%%%
% Run this code to convert lon coords if they use E/W, 
% instead of absolute value system
metlon = metlon - 360;
%%%

% Can use a DEM to define the mask area
[mask1, R1] = geotiffread('/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_mask.tif');
R1mat = georefobj2mat(R1, 'LL');
[masklon1, masklat1] = pixcenters(R1mat, size(mask1));
[masklon, masklat, ~] = grid2xyz(masklon1', masklat1', mask1);

% Use r.out.xyz to generate this from the basin mask raster
% maskxyz = dlmread('/Volumes/HD3/SWOTDA/Data/UMRB_2018/Basin/basincoords.txt', '|');
% masklon = maskxyz(maskxyz(:,3) == 1,1);
% masklat = maskxyz(maskxyz(:,3) == 1,2);
ncells = length(masklon);

% masklat = list of latitudes of the nonzero mask pixels. Must have four
% decimal places.

% kk = find(abs(masklat - 47.7813) < 0.0001); % these are the problematic pixels

% figure(3), imagesc(masklon, masklat, mask)

lat_ind = zeros(ncells,1);
lon_ind = zeros(ncells,1);
for k=1:ncells
        [~, lat_ind(k)] = min(abs(masklat(k) - metlat));
        [~, lon_ind(k)] = min(abs(masklon(k) - metlon));
end

%% Extract the forcings whose lat/lon match those of the basin mask

% Takes about 30 minutes to run for UMRB, five years.

disp('Extracting forcing variables')

nyears = endyear - beginyear + 1;
t_ind = 1;
cum_days = 0;

for t = beginyear:endyear
    
    if t==beginyear, tic, end
    prec = ncread(['prec.' num2str(t) '.nc'], 'prec');
    tmax = ncread(['tmax.' num2str(t) '.nc'], 'tmax');
    tmin = ncread(['tmin.' num2str(t) '.nc'], 'tmin');
    wind = ncread(['wind.' num2str(t) '.nc'], 'wind');
       
    info = ncinfo(['prec.' num2str(t) '.nc']);
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells);
        
    for k=1:ncells
        % Get the index of the met. forcing data that matches the basin mask
%         [Lia,lat_ind] = ismember(masklat(k),metlat);
%         [Lia,lon_ind] = ismember(masklon(k),metlon);
               
        data(:,1,k) = prec(lon_ind(k),lat_ind(k), :);        
        data(:,2,k) = tmin(lon_ind(k),lat_ind(k), :);
        data(:,3,k) = tmax(lon_ind(k),lat_ind(k), :);
        data(:,4,k) = wind(lon_ind(k),lat_ind(k), :);
    end
    
    if t_ind~=1
        data_cum = vertcat(data_cum, data);
    else
        data_cum = data;
    end

    t_ind = t_ind + 1;
    cum_days = size(prec,3) + cum_days;
    
    if t==beginyear 
        disp(['About ' num2str(toc*nyears/60) ...
            ' minutes remaining.'])
    end
    
    disp(t)
    
end

%%
% % Save met. forcings as .mat file
% % The dimensions are [numdays, numforcings, ncells]
% save('METFORC.mat', 'data_cum');

fstring = ['%.' num2str(grid_decimal) 'f'];
for k=1:ncells     
%     filename = ['data_' num2str(masklat(k),fstring) '_' num2str(masklon(k),fstring)];
    filename = ['data_' num2str(metlat(lat_ind(k)),fstring) '_' num2str(metlon(lon_ind(k)),fstring)];
    dlmwrite(fullfile(forcingsavedir, filename), data_cum(:,:,k), ' ')
end
display(['Forcing data saved to ' forcingsavedir])

%% Extract soils data

% Load the soils data

soils = load(fullfile(soilpath, soilname));

slat = soils(:,3);
slon = soils(:,4);

soils(:,1) = 0;

for k=1:ncells
    % Get the index of the met. forcing data that (nearly) matches the basin mask
%     sind = find(metlat(lat_ind(k)) == slat & metlon(lon_ind(k)) == slon);
        
%         [~, sind1] = min(abs(metlat(lat_ind(k)) - slat)); % gives the best latitude match
%         [~, sind2] = min(abs(metlon(lon_ind(k)) - slon)); % gives the best longitude match
        % need to find the best lat and lon match
        
        sind = find((abs(metlat(lat_ind(k)) - slat) <= 1/32) & (abs(metlon(lon_ind(k)) - slon) <= 1/32));
        % if the resolution is 1/16 degree, then we can use thresholds.
    soils(sind,1) = 1;
    if mod(k, 1000) == 0
        disp(k)
    end
end

soils1 = soils; % make a copy in case the next step goes sour
soils(soils(:,1) == 0,:) = []; % include this line to reduce file size

fstring = ['%.' num2str(grid_decimal) 'f'];
% fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];
% my 1/16 degree soil parameter file (global) with 2 soil layers
fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %d %d %.3f %.3f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d\n'];
fID = fopen(fullfile(soilsavedir, 'soils.global.umrb'),'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved to ' soilsavedir])

% Note: the delimiter and the format spec must be specified precisely as
% the Stehekin example from the VIC website in order to avoid the error
% about CELL LATITUDE not being found/for VIC to successfully read the soil
% parameter file.