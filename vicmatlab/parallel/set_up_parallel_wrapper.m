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

addpath(genpath('/home/jschap/Documents/Codes/VICMATLAB/vicmatlab'))

% /u/home/j/jschaper/vic/VICGlobal/v1_5/Classic/veglib_nh_new_resist.txt

%% 

% cd '/Volumes/HD3/SWOTDA'
cd ~/Documents/ESSD

% template soil parameter filename 
% control_params.soil_param = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_L13.txt';
% control_params.soil_param = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_soils_vg.txt';
control_params.soil_param = '/home/jschap/Documents/ESSD/clipped_soils_VG.txt';

% template global parameter filename
control_params.global_param = '/home/jschap/Documents/hoffman2transfer/global_param_WY2001-2011.txt';
% control_params.global_param = '/Volumes/HD4/SWOTDA/Data/UpperMiss/remote_par_classic_vg_1980-2011.txt';

% location where parallelized inputs for the VIC model should be saved
control_params.out_dir = '/home/jschap/Documents/hoffman2transfer/parallel/';
% control_params.out_dir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/parallel';

% full path for the VIC executable
control_params.vic_command = '/u/home/j/jschaper/vic/VIC-VIC.5.1.0.rc1/vic/drivers/classic/vic_classic.exe';

% location on Hoffman2 to write VIC outputs
control_params.vic_out_dir = '/u/scratch/j/jschaper/UCRB/WB_WY2001-2011_VG/';

% location on Hoffman2 from which global parameter files and soil parameter files will be read
% control_params.vic_in_dir = '/u/home/j/jschaper/vic/UMRB/EB_1980-2011/parallel';
control_params.vic_in_dir = '/u/scratch/j/jschaper/UCRB/WB_WY2001-2011_VG/parallel';

control_params.n_proc = 100; % number of processors
control_params.multidir = 0; % flag for using one or multiple directories for VIC outputs
control_params.soil_format = '3l'; % 3l or livneh

% cd '/Users/jschap/Documents/Codes/VICMATLAB/Control'
set_up_parallel(control_params)
% cd '/Volumes/HD3/SWOTDA'

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