% Subsets a Livneh et al. (2015) output variable to a particular time and
% place
%
% Written 1/29/2020 JRS
% Based on subset_livneh_runoff_baseflow.m
%
% Inputs
% addpath('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting/Subsetting')
% start_date = datetime(1979, 10, 1); % must start on the first day of the month
% end_date = datetime(2011, 9, 30);
% basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% fluxdir = '/Volumes/HD3/Livneh_2015/Fluxes';
% outname = '/Volumes/HD4/SWOTDA/Data/Colorado/L15';
% varname = 'SWE'
%
% Very much in development. Also a big time suck.

function [var_sub, timevector, masklon, masklat] = subset_livneh(start_date, end_date, basin_mask, fluxdir, outname, varname)

[mask1, R1, lon1, lat1] = geotiffread2(basin_mask);

% find the file containing the start date and read it
% continue reading files until you find the end date
% append data as you go

% Initialize arrays
ndays = days(end_date - start_date);
ncells = sum(sum(mask1>0));
var_sub = NaN(ndays, ncells);

rn_sub = NaN(ndays, ncells);

current_date = start_date;
first_iter = 1; % flag for first iteration
nt_old = 0;
d1 = 1;

%% First iteration to set up variables that will be used more than once

current_year = year(current_date);
current_month = month(current_date);

% current_month = 1;
if current_month < 10
    month_str = ['0' num2str(current_month)]; 
else
    month_str = num2str(current_month);
end
filename = ['Fluxes_Livneh_NAmerExt_15Oct2014.' num2str(current_year) month_str '.nc'];
filename = fullfile(fluxdir, filename);

var1 = ncread(filename, varname);
var1 = permute(var1, [2,1,3]);

lon_full = ncread(filename, 'lon');
lat_full = ncread(filename, 'lat');

xcoord = [lon1(1), lon1(end)];
ycoord = [lat1(1), lat1(end)];
figure, plotraster(xcoord, ycoord, mask1, 'Mask','Lon','Lat')

% This would be way easier to do in R. OK. Doing it in R.
[swe_cropped, croplon, croplat] = subset_netcdf_w_geotiffmask(var1(:,:,1), lon_full, lat_full, lon1, lat1, mask1);
swe_map = xyz2grid(masklon, masklat, swe_clipped');
figure, plotraster(xcoord, ycoord, swe_map, 'April SWE','Lon','Lat')


% figure, plotraster(xcoord, ycoord, var1(:,:,1), 'Oct. 1 SWE', 'Lon', 'Lat')

% Using net radiation to debug because it is easy to work with
rn = ncread(filename, 'NetRad');
rn = permute(rn, [2,1,3]);
    


figure
subplot(1,2,1)
plotraster(lon_full, lat_full, mean(rn,3), ['Net radiation (', datestr(current_date), ')'], 'Lon', 'Lat')
subplot(1,2,2)
plotraster(lon_full, lat_full, mean(var1,3), ['Var1 (', datestr(current_date), ')'], 'Lon', 'Lat')
caxis([0,1])

% get list of x, y, z coordinates for the basin mask
[mask1, R1, lon1, lat1] = geotiffread2(basin_mask);
[x,y,z] = grid2xyz(lon1', lat1', mask1); 
% put it back together (check)
mask2 = xyz2grid(x, y, z);
figure, plotraster(lon1, lat1, mask2, 'Reconstructed mask', '', '')

%  get list of x, y, z coordinates for the NetCDF data from Livneh
[xl, yl, zl] = grid2xyz(lon_full, lat_full, rn(:,:,1));
% put it back together (check)
rn2 = xyz2grid(xl, yl, zl);
figure, plotraster(lon_full, lat_full, rn2, 'Reconstructed Rn', '', '')

% Get the corresponding values of rn or of var1
B = [xl,yl];
[~, loc] = ismember([x, y],B,'rows');

% Some grid cells in the mask may not have exact matches in the Livneh domain
nodata = find(loc==0);
x(nodata) = [];
y(nodata) = [];
loc(loc==0) = [];

zl_sub = zl(loc);
rn3 = xyz2grid(x,y,zl_sub);

figure, plotraster(x, y, rn3, 'Reconstructed Rn', '', '')

%% Loop through all dates

while current_date < end_date
    
    current_year = year(current_date);
    current_month = month(current_date);
        
    if current_month < 10
        month_str = ['0' num2str(current_month)]; 
    else
        month_str = num2str(current_month);
    end
    filename = ['Fluxes_Livneh_NAmerExt_15Oct2014.' num2str(current_year) month_str '.nc'];
    filename = fullfile(fluxdir, filename);
    
    var1 = ncread(filename, varname);
    var1 = permute(var1, [2,1,3]);
    
    % Using net radiation to debug because it is easy to work with
    rn = ncread(filename, 'NetRad');
    rn = permute(rn, [2,1,3]);
            
    % Keep only the portion of the output that overlaps the study domain
    % Same methodology as vicinputworkflow for the subsetting    
        
    % Indices for days
    nt = size(var1, 3);
    d1 = d1 + nt_old;
    d2 = d1 + nt - 1;
    nt_old = nt;

    [rn_sub(d1:d2,:), masklon, masklat] = subset_netcdf_w_geotiffmask(rn, lon_full, lat_full, basin_mask);
    
    rn_map = xyz2grid(masklon, masklat, mean(rn_sub(d1:d2,:))');
    figure
    plotraster(masklon, masklat, rn_map, 'Net radiation', '', '')
    
    [var_sub(d1:d2,:), masklon, masklat] = subset_netcdf_w_geotiffmask(var1, lon_full, lat_full, basin_mask);
        
    current_date = current_date + nt;
    
    disp(current_date)
    
end

var_sub = var_sub';
timevector = start_date:end_date;
save(outname, 'var_sub', 'timevector');
disp(['Saved ' outname])

return