function OUTPUTS = check_outputs_wrapper(vic_out_dir)

% Performs various tasks to process the VIC Classic Driver model outputs
%
% Revised 10/7/2019 JRS
% Revised 3/11/2020 JRS, made into a function
%
% Note: if you do not want to use a parallel pool (e.g. if you don't have
% the license, then be sure to disable the automatic creation of parallel
% pools under your MATLAB preferences.

%% Set up the computing environment

timestep_out = 'daily';

wb_out_dir = fullfile(vic_out_dir, 'wb');
eb_out_dir = fullfile(vic_out_dir, 'eb');

% check that the VIC output files have been organized properly into /wb/ 
% and /eb/ folders
if exist(wb_out_dir, 'dir') == 0
    error('Make sure the VIC output files have been organized properly')
end

results_dir = fullfile(vic_out_dir, 'processed');
figdir = fullfile(vic_out_dir, 'figures');
mkdir(results_dir)
disp(['Created directory for results: ' results_dir]);
mkdir(figdir)
disp(['Created directory for figures: ' figdir]);

info = get_vic_run_metadata(wb_out_dir, eb_out_dir, timestep_out);

save(fullfile(results_dir, 'vic_run_metadata.mat'), 'info');
disp(['Saved VIC run metadata as ' fullfile(results_dir, 'vic_run_metadata.mat')])

%% Calculate time-average maps and area-average time series

% Using datastores to efficiently read in a large amount of data without
% needing to keep it all in memory at once. 
%
% However, this cannot be done in parallel, at least not without major 
% code adjustments, because using parfor here messes up the order of the
% pixels and gets the maps all jumbled up.

[wb_avg_map, wb_avg_ts, wb_sum_ts, ~, ~] = readVIC_ds(info.wb_out_dir, length(info.wbvars), info.ncells, info.nt);
[eb_avg_map, eb_avg_ts, eb_sum_ts, ~, ~] = readVIC_ds(info.eb_out_dir, length(info.ebvars), info.ncells, info.nt);

%% Assemble the data into an easy to deal with format

OUTPUTS = make_outputs_struct(info, wb_avg_ts, wb_avg_map, eb_avg_ts, eb_avg_map);
save(fullfile(results_dir, 'vic_outputs_summarized_daily.mat'), 'OUTPUTS');
disp(['Saved processed VIC outputs as ' fullfile(results_dir, 'vic_outputs_summarized_daily.mat')])

%% Make plots

% How to make a map:
% prec = wb_avg_map(:,13);
% prec_map = xyz2grid(info.lon, info.lat, prec);
% figure, plotraster(info.lon, info.lat, prec_map, 'precip', 'Lon', 'Lat', 1, gca)

plot_spatial_avg_ts(OUTPUTS, figdir);
disp(['Saved spatial average time series plots to ' figdir])

plot_time_avg_maps(OUTPUTS, figdir);
disp(['Saved time-average maps to ' figdir])

% %%

% OUTPUTS = make_outputs_struct(info, wb_avg_ts, wb_avg_map, eb_avg_ts, eb_avg_map);
% plot_time_avg_maps(OUTPUTS, figdir);

%% Write out GeoTiffs

% disp('Saving simulation-average water balance variables as geotiffs')
% xres = 1/16;
% yres = 1/16;
% 
% R = makerefmat(min(info.lon), min(info.lat), xres, yres);
% geotiffwrite(fullfile(figdir, 'average_temperature.tif'), OUTPUTS.avg_map.OUT_AIR_TEMP, R)
% geotiffwrite(fullfile(figdir, 'average_precipitation.tif'), OUTPUTS.avg_map.OUT_PREC, R)
% geotiffwrite(fullfile(figdir, 'average_evaporation.tif'), OUTPUTS.avg_map.OUT_EVAP, R)
% geotiffwrite(fullfile(figdir, 'average_runoff.tif'), OUTPUTS.avg_map.OUT_RUNOFF, R)
% geotiffwrite(fullfile(figdir, 'average_baseflow.tif'), OUTPUTS.avg_map.OUT_BASEFLOW, R)
% geotiffwrite(fullfile(figdir, 'average_swe.tif'), OUTPUTS.avg_map_snow.OUT_SWE, R)
% disp(['Wrote out Geotiffs of select water balance variables to ' figdir])

return
