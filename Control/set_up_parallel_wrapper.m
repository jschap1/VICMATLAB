% Parallel wrapper

addpath('/Volumes/HD3/SWOTDA/Calibration/CalVIC')

control_params = struct();
control_params.vic_out_dir = '';
control_params.soil_param = '';
control_params.global_param_file = '';
control_params.n_proc = 40;
control_params.vic_command = '';

set_up_parallel()