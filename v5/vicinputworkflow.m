function [] = vicinputworkflow()

% INPUTS
%
% Basin mask, as list of coordinates.
% Run the R script routinputworkflow first to generate the basin mask
% Use r.out.xyz to the list of coordinates generate this from the basin mask raster
%
% Designed for VIC 5, specifically with MERRA-2 forcings in mind
%
% Met. forcings at daily resolution, saved as separate NetCDF files for
% each variable and year (e.g. prec.1965.nc, wind.1990.nc)
% Basin mask with using appropriate resolution, geographic coordinates
% Soil parameter file covering >= the domain of interest
% Designed for use with Livneh met. forcing data. 
% Adapted from VIC 4 version 2/22/2019 JRS
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


%% Specify inputs

% cd '/Volumes/HD3/SWOTDA/Data/IRB/'

% Directory of daily CONUS met. forcing files
forcingpath = '/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Processed/MERRA2_Outputs';

% Name of basin mask at VIC modeling resolution
% basinmask = './Data/UMRB/Basin/basin.wgs.asc';
% basinmask = './Data/UMRB/rout/umrb.fract';

% Directory where clipped forcing files should be saved
forcingsavedir = './VIC/FORC';

% Number of forcings in the forcing file
numforcings = 7;

% The MERRA forcing files are hourly, and there is one set of forcings per
% year

% Beginning and ending years of simulation (must be included in the daily CONUS
% met. forcing file)
beginyear = 1990;
endyear = 2000;

% Number of decimal points of precision to use for forcing file names
grid_decimal = 5;

% Directory of CONUS soil parameter file
soilpath = '/Volumes/HD3/SWOTDA/Data/Global/second_attempt';
soilname = 'global_soils_1_16.txt'; 

% Directory where clipped soil parameter file should be saved
soilsavedir = './VIC';

%% Load the mask

% addpath(forcingpath)
% metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
% metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');

metlat = 1; % list of coordinates where forcing data are available
metlon = 1; 
% Ask GK to get this for me...

fnames = dir(fullfile(forcingpath, '*.out'));

% Get list of pixel numbers --------------------
n_merra = length(fnames);
pixelno = zeros(n_merra,1);
for k=1:n_merra
    temp1 = split(fnames(k).name, '_year');
    temp2 = split(temp1{1}, 'pixel_');
    pixelno(k) = str2double(temp2{2});
end

% Use static file to get lat and lon from pixel numbers
static_file = load('/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Processed/static_file_MERRA2.in');
% [~, basin_inds] = ismember(pixelno, static_file(:,1));
metlat = static_file(pixelno,2);
metlon = static_file(pixelno,3);
% ----------------------------------------------

%%%
% Run this code to convert lon coords if they use E/W, 
% instead of absolute value system
% metlon = metlon - 360;
%%%

% Use r.out.xyz to generate this from the basin mask raster
maskxyz = dlmread('basincoords.txt', '|');
masklon = maskxyz(maskxyz(:,3) == 1,1);
masklat = maskxyz(maskxyz(:,3) == 1,2);
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
    
end

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
    % Get the index of the met. forcing data that matches the basin mask
    sind = find(metlat(lat_ind(k)) == slat & metlon(lon_ind(k)) == slon);
    soils(sind,1) = 1;
end

soils(soils(:,1) == 0,:) = []; % include this line to reduce file size

fstring = ['%.' num2str(grid_decimal) 'f'];
fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];
fID = fopen(fullfile(soilsavedir, 'soils.SB'),'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved to ' soilsavedir])

% Note: the delimiter and the format spec must be specified precisely as
% the Stehekin example from the VIC website in order to avoid the error
% about CELL LATITUDE not being found/for VIC to successfully read the soil
% parameter file.