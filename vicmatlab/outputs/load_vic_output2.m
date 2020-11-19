% Function for loading one particular VIC variable, as a time series
% 
% Inspired by subset_livneh.m function for image mode variables, but this
% version is for classic mode outputs. Code also draws heavily from
% VICMATLAB::check_outputs_wrapper.m
%
% Revised 8/4/2020 to be more flexible w/regard to inputs
%
% Inputs
% varargin = timestep_out, saveflag, savename

function [timevector, var1, info] = load_vic_output2(vic_out_dir, varname, prefix, varargin)

%% Parse optional arguments;

numvarargs = length(varargin);
if numvarargs > 4
    error('The max number of optional arguments is 4')
end

optargs = {3, 'daily', 0, './varname.mat'};
optargs(1:numvarargs) = varargin;
[headerlines, timestep_out, saveflag, savename] = optargs{:};

info = get_vic_run_metadata2(vic_out_dir, timestep_out, prefix, headerlines);

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

official_variable_name = lutable{2}{var_index};
var_col = find(strcmp(official_variable_name, info.vars));

%% Load VIC output variable

% Initialize arrays
varnames = dir([vic_out_dir '/' prefix '*.txt']);
dat = dlmread(fullfile(vic_out_dir, varnames(1).name), '\t', 3, 0);

timevector = datetime(dat(:,1), dat(:,2), dat(:,3));

ndays = size(dat, 1);
ncells = length(varnames);
problem_flag = 0;

var1 = zeros(ncells, ndays);
for k=1:ncells
    dat = dlmread(fullfile(vic_out_dir, varnames(k).name), '\t', 3, 0);
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

%% Save output (optional)

if saveflag
    save(savename, 'timevector', 'var1', 'info');
    disp(['Saved variable to ' savename])
end

return
