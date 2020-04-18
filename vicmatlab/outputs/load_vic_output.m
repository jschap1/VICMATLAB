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

function [timevector, var1, info] = load_vic_output(vic_out_dir, varname, varargin)
% function [timevector, var_sub, var1, info] = load_vic_output(vic_out_dir, basin_mask, var_col)

%% Parse optional arguments;

numvarargs = length(varargin);
if numvarargs > 3
    error('The max number of optional arguments is 3')
end

optargs = {'daily', 0, './varname.mat'};
optargs(1:numvarargs) = varargin;
[timestep_out, saveflag, savename] = optargs{:};

info = get_vic_run_metadata(vic_out_dir, timestep_out);

% [mask1, ~, lon1, lat1] = geotiffread2(basin_mask);
% figure, plotraster(lon1, lat1, mask1, '','','')

%% Parse the input variable name:

fID = fopen('vic_output_lookup_table');
lutable = textscan(fID, '%s%s%s');
fclose(fID);

% Search the table for a match
var_index = find(strcmp(varname, lutable{1}));
if isempty(var_index)
    var_index = find(strcmp(varname, lutable{2}));
end
if isempty(var_index)
    error('Variable not found. Check spelling.')
end

ebwb = lutable{3}{var_index};
official_variable_name = lutable{2}{var_index};
if strcmp(ebwb, 'wb')
    out_dir = fullfile(vic_out_dir, 'wb');
    var_col = find(strcmp(official_variable_name, info.wbvars));
else
    out_dir = fullfile(vic_out_dir, 'eb');
    var_col = find(strcmp(official_variable_name, info.ebvars));
end

%% Load VIC output variable

% Initialize arrays
varnames = dir([out_dir '/*.txt']);
dat = dlmread(fullfile(out_dir, varnames(1).name), '\t', 3, 0);

% size(dat) % ndays by 27 for VICGlobal case

timevector = datetime(dat(:,1), dat(:,2), dat(:,3));

ndays = size(dat, 1);
ncells = length(varnames);
problem_flag = 0;

var1 = zeros(ncells, ndays);
for k=1:ncells
    dat = dlmread(fullfile(out_dir, varnames(k).name), '\t', 3, 0);
    try
        var1(k,:) = dat(:,var_col);
    catch
        % Some of the grid cells may not have run to completion for
        % whatever reason. If this happens, check to see what the problem
        % is with the VIC output files.
        if problem_flag == 0
            problem_flag = 1;
            disp('Some of the grid cells may not have run to completion')
        end
        nn = size(dat(:,var_col));
        var1(k,1:nn) = dat(:,var_col);        
    end
    % 1 by ndays and ndays by 1 for VICGlobal case
    if mod(k,1000)==0
        disp(k) % update this to a proper progress bar
    end
end

%% Subset the variable (not implemented)

% % There may or may not be a need to subset var1.
% % 
% % It depends if the VIC simulation was run over the study basin, precisely,
% % or if the simulation was run over a larger area and we wish to subset the
% % results to the basin.
% % 
% % The following code checks which case is happening and performs the
% % appropriate computations:
% % 
% % basinmask;
% % 
% % mask_vect = mask1(:);
% % var_sub = var1;
% % 
% % if size(var_sub, 1) < length(mask_vect)
% %     disp('The number of VIC grid cells is smaller than or equal to the number of grid cells in the basin mask')
% % else
% %     disp('The number of VIC grid cells is larger than the number of grid cells in the basin mask')
% %     disp('The user should crop the output to the basin mask.')
% %     disp('Attempting to crop the results to the basin mask')
% %     var_sub = var_sub(mask_vect==1,:);
% %     var_sub(mask_vect ~= 1,:) = [];    
% % end
% % 
% % Use the mask to subset var1. See subset_netcdf_w_geotiff_mask.m from
% % subset_livneh.m. Ideally, there is no need to call xyz2grid, and we can
% % preserve the time-resolution of the data.
% % 
% % [var_sub, ~, ~] = subset_netcdf_w_geotiffmask(var1, lon1, lat1, basin_mask);
% % The above would work if var_sub were nx by ny by nt grid
% % 
% % Will need to find another way (1/29/2020)
% % 
% % ncells_in_mask = sum(sum(mask1>0));
% % var_sub = NaN(ndays, ncells_in_mask);

%% Save output (optional)

if saveflag
    save(savename, 'timevector', 'var1', 'info');
    disp(['Saved variable to ' savename])
end

return
