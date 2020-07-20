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
%
% Update 11/25/2019 JRS to be generalized, not just call-able from the
% CalVIC workflow, as it was previously. Also removed some extraneous code.
%
% INPUTS
% control_params -- a structure with the following fields:
%
% soil_param = template soil parameter filename 
% global_param = template global parameter filename
%
% out_dir = location where parallelized inputs for the VIC model should be saved
% vic_command = full path for the VIC executable
%
% vic_out_dir = location on Hoffman2 to write VIC outputs
% vic_in_dir = location on Hoffman2 from which global parameter files and soil parameter
% files will be read
%
% n_proc = number of processors
% multidir = flag for using one or multiple directories for VIC outputs
% soil_format = 3l or livneh

function outname = set_up_parallel(control_params)

%% Inputs

% Specify lines of the global parameter file for the following:
state_line_ind = 111; % 110, 111
log_line_ind = 182; % 181, 195
out_line_ind = 183; % 182, 196
soil_line_ind = 154; % 153, 167

mkdir(control_params.out_dir)
disp(['Created directory ' control_params.out_dir ' to store the subsetted soil parameter files'])

[~, soil_basename, ~] = fileparts(control_params.soil_param);
[~, global_basename, ~] = fileparts(control_params.global_param);

soils_all = dlmread(control_params.soil_param);
ncells = size(soils_all, 1);

n = control_params.n_proc;
nl = ceil(ncells/n); % number of grid cells per soil parameter file
disp(['Each soil parameter file will have ' num2str(nl) ' grid cells'])

% Split into n directories so the outputs from each division are saved in different directories (for debugging purposes)
if control_params.multidir
    disp('Creating separate directories for outputs')
    outdir_names = cell(n,1);
    for jj=1:n
      outdir_names{jj} = fullfile(control_params.vic_out_dir, ['Output' num2str(jj)]);
    end
else
    outdir_names = control_params.vic_out_dir;
end

%% Make copies of the global parameter file

disp('Copying the global parameter file')
for jj=1:n

    % read the global parameter file
    A = read_global_param_file(control_params.global_param);

    % change the soil parameter line
    
    A{soil_line_ind} = ['SOIL ' fullfile(control_params.vic_in_dir, [soil_basename, '_', num2str(jj), '.txt'])];
   
    % change the output directory
    if control_params.multidir
        A{log_line_ind} = ['LOG_DIR ' outdir_names{jj} '/'];
        A{out_line_ind} = ['RESULT_DIR ' outdir_names{jj} '/'];
%         A{state_line_ind} = ['STATENAME ' outdir_names{jj} '/state_file'];
    else
        A{log_line_ind} = ['LOG_DIR ' outdir_names '/'];
        A{out_line_ind} = ['RESULT_DIR ' outdir_names '/'];
%         A{state_line_ind} = ['STATENAME ' outdir_names '/state_file'];
    end

    sn = fullfile(control_params.out_dir, [global_basename, '_', num2str(jj), '.txt']);
    write_global_param_file(A, sn)

end

%% Subdivide the soil parameter file

disp('Sub-dividing the soil parameter file')

for jj=1:n
     
    if jj==n
        soils_sub = soils_all((jj-1)*nl+1:ncells,:);
    else
        soils_sub = soils_all((jj-1)*nl+1:(jj-1)*nl+nl,:);
    end
    
    outname = fullfile(control_params.out_dir, [soil_basename, '_', num2str(jj) '.txt']);
    precision = 5;
    write_soils(precision, soils_sub, outname, control_params.soil_format);

end

%% Make the VIC parallel wrapper

[globpath, globname globext] = fileparts(control_params.global_param);
% parfilenam = fullfile(globpath, [globname '_${SGE_TASK_ID}' globext]);
parfilenam = fullfile(control_params.vic_in_dir, [globname '_${SGE_TASK_ID}' globext]);
vic_run_command = [control_params.vic_command ' -g ' parfilenam];

outname = fullfile(fileparts(control_params.out_dir), 'vic_parallel_wrapper.sh');
A = cell(3,1);
A{1} = '#!/bin/bash';
A{2} = 'echo $SGE_TASK_ID';
A{3} = vic_run_command;
A{4} = -1;
write_global_param_file(A, outname)
disp(['Wrote exec file to ' outname])

return
