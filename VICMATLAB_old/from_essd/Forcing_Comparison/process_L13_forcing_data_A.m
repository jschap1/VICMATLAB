function forc = process_L13_forcing_data_A(forcingpath, start_date, end_date, target_lon, target_lat)

% Adapted from VICMATLAB/v4/vicinputworkflow.m to load forcing data for a
% particular grid cell from the NetCDF file
%
% INPUTS
% forcingpath = directory of forcing files (netCDF)
% start_date
% end_date
% target_lon
% target_lat
%
% OUTPUTS
% precipitation
% tmax
% tmin
% wind speed

%% Load the mask

% Directory of daily CONUS met. forcing file
% forcingpath = '/Volumes/HD3/Livneh_2013/MetNC';

% Beginning and ending years of simulation (must be included in the daily CONUS
% met. forcing file)
beginyear = year(start_date);
endyear = year(end_date);

addpath(forcingpath)
metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');

%%%
% Run this code to convert lon coords if they use E/W, 
% instead of absolute value system
metlon = metlon - 360;
%%%

% Can use a DEM to define the mask area
% [mask1, R1] = geotiffread('/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_mask.tif');
% R1mat = georefobj2mat(R1, 'LL');
% [masklon1, masklat1] = pixcenters(R1mat, size(mask1));
% [masklon, masklat, ~] = grid2xyz(masklon1', masklat1', mask1);

% Use r.out.xyz to generate this from the basin mask raster
% maskxyz = dlmread('/Volumes/HD3/SWOTDA/Data/UMRB_2018/Basin/basincoords.txt', '|');
% masklon = maskxyz(maskxyz(:,3) == 1,1);
% masklat = maskxyz(maskxyz(:,3) == 1,2);
% ncells = length(masklon);

% masklat = list of latitudes of the nonzero mask pixels. Must have four
% decimal places.

% kk = find(abs(masklat - 47.7813) < 0.0001); % these are the problematic pixels

% figure(3), imagesc(masklon, masklat, mask)

[~, lat_ind] = min(abs(target_lat - metlat));
[~, lon_ind] = min(abs(target_lon - metlon));

%% Extract the forcings whose lat/lon match

disp('Extracting forcing variables')

start1 = [lon_ind,lat_ind,1];
count1 = [1,1,Inf];

for t = beginyear:endyear
    
    if t==beginyear, tic, end
    
    prec = ncread(['prec.' num2str(t) '.nc'], 'prec', start1, count1);
    tmax = ncread(['tmax.' num2str(t) '.nc'], 'tmax', start1, count1);
    tmin = ncread(['tmin.' num2str(t) '.nc'], 'tmin', start1, count1);
    wind = ncread(['wind.' num2str(t) '.nc'], 'wind', start1, count1);
    
    if t==beginyear
        forc = struct();
        forc.prec = squeeze(prec);
        forc.tmax = squeeze(tmax);
        forc.tmin = squeeze(tmin);
        forc.wind = squeeze(wind);
    else
        forc.prec = vertcat(forc.prec, squeeze(prec));
        forc.tmax = vertcat(forc.tmax, squeeze(tmax));
        forc.tmin = vertcat(forc.tmin, squeeze(tmin));
        forc.wind = vertcat(forc.wind, squeeze(wind));        
    end
           
    disp(['Finished processing calendar year ' num2str(t)])
    
end

forc.time_daily = (datetime(beginyear, 1, 1):datetime(endyear, 12, 31))';

return
