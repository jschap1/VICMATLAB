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

function outname = set_up_parallel(control_params)

%% Inputs

% global_param_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/globalparam_WB_SB_1980-2018.txt';
% soil_parameter_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/soils_snapped.SB';

savedir = fullfile(control_params.vic_out_dir, 'parallel');
mkdir(savedir)
% savedir_soil = savedir; % directory to save the subsetted soil param files
% savedir = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/vic_parallel_WB_1980-2018'; % directory to save the newly created global param files

[~, soil_basename, ~] = fileparts(control_params.soil_param);
[~, global_basename, ~] = fileparts(control_params.global_param_file);

% global_basename = './global_param_irb';

soils_all = dlmread(control_params.soil_param);
ncells = size(soils_all, 1);

% how many divisions to divide the soil parameter file into
n = control_params.n_proc;
nl = ceil(ncells/n); 

% ncells = 75231;
% nl = 1000; % number of grid cells to keep in each soil parameter file
% n = ceil(ncells/nl); 

% Create n directories so the outputs from each division are saved in different directories (for debugging purposes)
% outdir_names = cell(n,1);
% for jj=1:n
%   outdir_names{jj} = fullfile(savedir_out, ['Output_' num2str(jj)]);
% %   mkdir(outdir_names{jj})
% end

%% Make copies of the global parameter file

for jj=1:n

    % read the global parameter file
    fid = fopen(control_params.global_param_file, 'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);

    % change the soil parameter line
    A{153} = ['SOIL ' fullfile(savedir, [soil_basename, '_', num2str(jj), '.txt'])];

    % change the output directory
    A{181} = ['LOG_DIR ' control_params.vic_out_dir];
    A{182} = ['RESULT_DIR ' control_params.vic_out_dir];

    % Write cell A into txt
    sn = fullfile(savedir, [global_basename, '_', num2str(jj), '.txt']);
    
    fid = fopen(sn, 'w');
    for i=1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

end

%% Subdivide the soil parameter file

for jj=1:n
     
    if jj==n
        soils_sub = soils_all((jj-1)*nl+1:ncells,:);
    else
        soils_sub = soils_all((jj-1)*nl+1:(jj-1)*nl+nl,:);
    end
    
    outname = fullfile(savedir, [soil_basename, '_', num2str(jj) '.txt']);
    precision = 5;
    write_soils(precision, soils_sub, outname, '3l');

end

%% Make the VIC parallel wrapper

[globpath, globname globext] = fileparts(control_params.global_param_file);
% parfilenam = fullfile(globpath, [globname '_${SGE_TASK_ID}' globext]);
parfilenam = fullfile(savedir, [globname '_${SGE_TASK_ID}' globext]);
vic_run_command = [control_params.vic_command ' -g ' parfilenam];

outname = fullfile(fileparts(control_params.vic_out_dir), 'vic_parallel_wrapper.sh');
A = cell(3,1);
A{1} = '#!/bin/bash';
A{2} = 'echo $SGE_TASK_ID';
A{3} = vic_run_command;
A{4} = -1;
write_global_param_file(A, outname)
disp(['Wrote exec file to ' outname])

return
