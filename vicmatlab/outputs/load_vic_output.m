% Function for loading one particular VIC variable, as a time series
% 
% Can only read data from the wb output file, as written.
%
% Inspired by subset_livneh.m function for image mode variables, but this
% version is for classic mode outputs. Code also draws heavily from
% VICMATLAB::check_outputs_wrapper.m
%
% Inputs
% var_col = column of the variable you wish to extract (use 27 for SWE)

function [timevector, var_sub, var1, info] = load_vic_output(vic_out_dir, basin_mask, var_col)

timestep_out = 'daily';

wb_out_dir = fullfile(vic_out_dir, 'wb');
eb_out_dir = fullfile(vic_out_dir, 'eb');

info = get_vic_run_metadata(vic_out_dir, timestep_out);

[mask1, ~, lon1, lat1] = geotiffread2(basin_mask);
% figure, plotraster(lon1, lat1, mask1, '','','')

% Initialize arrays
wbnames = dir([wb_out_dir '/*.txt']);
dat = dlmread(fullfile(wb_out_dir, wbnames(1).name), '\t', 3, 0);
timevector = datetime(dat(:,1), dat(:,2), dat(:,3));

ndays = size(dat, 1);
ncells = length(wbnames);

var1 = zeros(ncells, ndays);
for k=1:ncells
    dat = dlmread(fullfile(wb_out_dir, wbnames(k).name), '\t', 3, 0);
    var1(k,:) = dat(:,var_col);
    if mod(k,1000)==0
        disp(k)
    end
end

% There may or may not be a need to subset var1.

% It depends if the VIC simulation was run over the study basin, precisely,
% or if the simulation was run over a larger area and we wish to subset the
% results to the basin.

% The following code checks which case is happening and performs the
% appropriate computations:

mask_vect = mask1(:);
var_sub = var1;

if size(var_sub, 1) < length(mask_vect)
    disp('The number of VIC grid cells is smaller than or equal to the number of grid cells in the basin mask')
else
    disp('The number of VIC grid cells is larger than the number of grid cells in the basin mask')
    disp('Attempting to crop the results to the basin mask')
    var_sub = var_sub(mask_vect==1,:);
    var_sub(mask_vect ~= 1,:) = [];    
end

% Use the mask to subset var1. See subset_netcdf_w_geotiff_mask.m from
% subset_livneh.m. Ideally, there is no need to call xyz2grid, and we can
% preserve the time-resolution of the data.

% [var_sub, ~, ~] = subset_netcdf_w_geotiffmask(var1, lon1, lat1, basin_mask);
% The above would work if var_sub were nx by ny by nt grid
%
% Will need to find another way (1/29/2020)

% ncells_in_mask = sum(sum(mask1>0));
% var_sub = NaN(ndays, ncells_in_mask);

return
