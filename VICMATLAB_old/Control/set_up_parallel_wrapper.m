% Parallel wrapper
%
% Written 11/25/2019 JRS
%
% Sets up the VIC global parameter and soil parameter files for a
% multi-processor run on Hoffman2. 
%
% Eventually, it would be cool to make a script that completely prepares
% the global parameter file from scratch based on a Matlab wrapper,
% potentially with a GUI, like Sungwook Wi did in VIC-ASSIST.

addpath('/Users/jschap/Documents/Codes/VICMATLAB/Control')
addpath('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting')
addpath('/Users/jschap/Documents/Codes/VICMATLAB/Make_Soils')

%% 

cd '/Volumes/HD3/SWOTDA'

% template soil parameter filename 
control_params.soil_param = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_L13.txt';
% control_params.soil_param = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_soils_vg.txt';

% template global parameter filename
control_params.global_param = '/Volumes/HD4/SWOTDA/Data/Colorado/remote_par_classic_l13_2000-2011.txt';
% control_params.global_param = '/Volumes/HD4/SWOTDA/Data/UpperMiss/remote_par_classic_vg_1980-2011.txt';

% location where parallelized inputs for the VIC model should be saved
control_params.out_dir = '/Volumes/HD4/SWOTDA/Data/Colorado/parallel_L13';
% control_params.out_dir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/parallel';

% full path for the VIC executable
control_params.vic_command = '/u/home/j/jschaper/vic/VIC-VIC.5.1.0.rc1/vic/drivers/classic/vic_classic.exe';

% location on Hoffman2 to write VIC outputs
control_params.vic_out_dir = '/u/scratch/j/jschaper/UCRB/EB_2000-2011_L13/';

% location on Hoffman2 from which global parameter files and soil parameter files will be read
% control_params.vic_in_dir = '/u/home/j/jschaper/vic/UMRB/EB_1980-2011/parallel';
control_params.vic_in_dir = '/u/home/j/jschaper/vic/UCRB/EB_2000-2011_L13/parallel';

control_params.n_proc = 100; % number of processors
control_params.multidir = 1; % flag for using one or multiple directories for VIC outputs
control_params.soil_format = 'livneh'; % 3l or livneh

cd '/Users/jschap/Documents/Codes/VICMATLAB/Control'
set_up_parallel(control_params)
cd '/Volumes/HD3/SWOTDA'

% Remember to create output directories on Hoffman2 before executing the
% script. 

%% Sample inputs

% control_params = struct();
% % control_params.vic_out_dir = '/Volumes/HD4/SWOTDA/Data/IRB/EB_JF1980';
% control_params.vic_out_dir = '/u/scratch/j/jschaper/EB_JF1980';
% control_params.soil_param = '/Volumes/HD4/SWOTDA/Data/IRB/soils_irb.txt';
% control_params.global_param_file = '/Volumes/HD4/SWOTDA/Data/IRB/par_EB_JF1980.txt';
% control_params.n_proc = 40;
% control_params.vic_command = '/u/home/j/jschaper/vic/VIC-VIC.5.1.0.rc1/vic/drivers/classic/vic_classic.exe';

% control_params = struct();
% control_params.vic_out_dir = '/u/scratch/j/jschaper/EB_JF1980';
% control_params.soil_param = '/u/home/j/jschaper/vic/IRB/irb_soils.txt';
% control_params.global_param_file = '/u/home/j/jschaper/vic/IRB/par_EB_JF1980.txt';
% control_params.n_proc = 40;
% control_params.vic_command = '';