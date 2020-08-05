% Set up parallel
%
% Subsets the VIC soil parameter file into multiple, smaller soil parameter files 
% and generates corresponding copies of the VIC global parameter file for running 
% VIC 4 or VIC 5 Classic in parallel as a job array on Hoffman2. 
% This works because computations for one grid cell are independent of other grid cells.
%
% Written by Dongyue Li 4/2019
% Modified 5/1/2019 JRS
% Updated 11/18/2019 JRS to be a function, not a script
% Updated 7/18/2020 JRS to use a simpler directory structure
%
% INPUTS
% outdir = Directory where VIC should save outputs (the parallel inputs will be
% written here, too)
% soil_param = Name of soil parameter file
% global_param_file = Name of global parameter file
% n_proc = Number of processors to use
% vic_command = command to call VIC

function outname = set_up_parallel(control_params)

%% Inputs

mkdir(control_params.vic_outdir);
mkdir(control_params.log_outdir);
mkdir(control_params.vic_indir);

[~, soil_basename, ~] = fileparts(control_params.soil_param);
[~, global_basename, ~] = fileparts(control_params.global_param_file);

soils_all = dlmread(control_params.soil_param);
ncells = size(soils_all, 1);

% how many divisions to divide the soil parameter file into
n = control_params.n_proc;
nl = ceil(ncells/n); 

%% Make copies of the global parameter file

% Lines in the global parameter file where the soil parameters, results,
% and log inputs are specified. Detect automatically.

soil_line = find_line(control_params.global_param_file, 'Soil parameter path/file');
result_line = find_line(control_params.global_param_file, 'RESULT_DIR');
log_line = find_line(control_params.global_param_file, 'LOG_DIR');

for jj=1:n

    A = read_global_param_file(control_params.global_param_file);
    
    % change the soil parameter line
    A{soil_line} = ['SOIL ' fullfile(control_params.vic_indir, [soil_basename, '_', num2str(jj), '.txt'])];

    % change the output directory
    A{log_line} = ['LOG_DIR ' control_params.log_outdir];
    A{result_line} = ['RESULT_DIR ' control_params.vic_outdir];

    % Write out the modified global parameter file
    sn = fullfile(control_params.vic_indir, [global_basename, '_', num2str(jj), '.txt']);
    write_global_param_file(A, sn)

end

%% Subdivide the soil parameter file

for jj=1:n
     
    if jj==n
        soils_sub = soils_all((jj-1)*nl+1:ncells,:);
    else
        soils_sub = soils_all((jj-1)*nl+1:(jj-1)*nl+nl,:);
    end
    
    outname = fullfile(control_params.vic_indir, [soil_basename, '_', num2str(jj) '.txt']);
    precision = 5;
    write_soils(precision, soils_sub, outname, control_params.soil_format);

end

%% Make the VIC parallel wrapper

[globpath, globname globext] = fileparts(control_params.global_param_file);
% parfilenam = fullfile(globpath, [globname '_${SGE_TASK_ID}' globext]);
parfilenam = fullfile(control_params.vic_indir, [globname '_${SGE_TASK_ID}' globext]);
vic_run_command = [control_params.vic_command ' -g ' parfilenam];

outname = fullfile(fileparts(control_params.vic_indir), 'vic_parallel_wrapper.sh');
A = cell(3,1);

% Different command file for local computer vs. Hoffman2
if control_params.local
    A{1} = '#!/bin/bash';
    A{2} = ['for SGE_TASK_ID in {1..', num2str(control_params.n_proc), '}'];
    A{3} = 'do';
    A{4} = 'echo $SGE_TASK_ID';
    A{5} = [vic_run_command, ' &']; % run in background
    A{6} = 'done';
    A{7} = -1;
else
    A{1} = '#!/bin/bash';
    A{2} = 'echo $SGE_TASK_ID';
    A{3} = vic_run_command;
    A{4} = -1;
end

write_global_param_file(A, outname)
system(['chmod +x ' outname])
disp(['Wrote exec file to ' outname])

return

%% Scrap 

% global_param_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/globalparam_WB_SB_1980-2018.txt';
% soil_parameter_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/soils_snapped.SB';
% savedir_soil = savedir; % directory to save the subsetted soil param files
% savedir = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/vic_parallel_WB_1980-2018'; % directory to save the newly created global param files
% global_basename = './global_param_irb';
% ncells = 75231;
% nl = 1000; % number of grid cells to keep in each soil parameter file
% n = ceil(ncells/nl); 
% Create n directories so the outputs from each division are saved in different directories (for debugging purposes)
% outdir_names = cell(n,1);
% for jj=1:n
%   outdir_names{jj} = fullfile(savedir_out, ['Output_' num2str(jj)]);
% %   mkdir(outdir_names{jj})
% end