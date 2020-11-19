% Load Classic Mode Inputs
%
% 8/5/2020 JRS
% Loads classic mode outputs from VIC-5 Classic
% 
% INPUTS
% vic_out_dir = directory where raw VIC outputs are stored
% prefix = prefix of VIC output files (e.g. 'fluxes'). Must include '_'
% timestep_out = either 'daily' or 
% processed_dir = directory to store processed VIC outputs
% basinmask = name of basin mask (geotiff)
%
% OUTPUT
%
% classic = structure containing all data from VIC simulation
% This variable will be very large for a big VIC simulation
%
% For simulations with large numbers of grid cells and/or time steps, avoid
% using load_classic_mode_outputs. Instead, refer to ucrb_wrapper.m or
% similar. Can use load_outputs as toggle here. See below.

function classic = load_classic_mode_outputs(vic_out_dir, prefix, timestep_out, processed_dir, basinmask)

mkdir(processed_dir); % only need to do this if it doesn't already exist.

% Post spin-up, original parameter simulation
classic.outdir = vic_out_dir;

classic.metadata = get_vic_run_metadata2(classic.outdir, timestep_out, prefix);
classic.processed_dir = processed_dir;
[basin_mask, Rmask, masklon, masklat] = geotiffread2(basinmask);

if ~strcmp(class(basin_mask), 'double')
    basin_mask = double(basin_mask);
end

% starting index at 4 because the first three are year, month, and day
for i=4:length(classic.metadata.vars)
    varname = classic.metadata.vars{i};
    savename = fullfile(classic.processed_dir, [varname '.mat']);
    try
        [classic.timevector, ~, classic.info] = load_vic_output2(classic.outdir, varname, prefix, 'daily', 1, savename);
        masksavename = fullfile(classic.processed_dir, [varname '_masked.mat']);
        mask_results_to_bb(savename, basin_mask, masksavename);
    catch
        disp(['Failed to process ' varname])
        disp(savename)
        disp(masksavename)
    end
end
 
load_outputs = 1; % toggle for loading processed outputs into MATLAB

% Load classic outputs
if load_outputs
    for i=4:length(classic.metadata.vars)
        varname = classic.metadata.vars{i};
        tempname = strsplit(varname, 'OUT_');
        masksavename = fullfile(classic.processed_dir, [varname '_masked.mat']);
        temp = load(masksavename);
        classic.out.(tempname{2}) = temp.output_map_masked;
    end
end

return