% Get VIC run metadata
%
% Gets metadata from the VIC run, such as simulation start and end times,
% names of output variables, etc.
%
% Inputs
% vic_out_dir
% prefix
% timestep_out = 'daily' or 'hourly'
%
% Adapted from get_vic_run_metadata.m to allow different output filenames
% besides 'eb' and 'wb'

function info = get_vic_run_metadata2(vic_out_dir, timestep_out, prefix)

headerlines = 3;

% Get metadata (lat, lon, time, names)
fluxnames = dir(fullfile(vic_out_dir, [prefix, '*']));

[lon, lat] = get_coordinates_from_VIC_file(vic_out_dir, prefix);

sample_flux_file = fullfile(vic_out_dir, fluxnames(1).name);
flux_varnames = get_vic_header(sample_flux_file, headerlines);

flux_out = dlmread(sample_flux_file, '\t', headerlines, 0);
switch timestep_out
    case 'hourly'
        date_array = flux_out(:,1:4);
        timevector = datetime(date_array(:,1), date_array(:,2), date_array(:,3), 0, 0, date_array(:,4));
    case 'daily'
        date_array = flux_out(:,1:3);
        timevector = datetime(date_array(:,1), date_array(:,2), date_array(:,3));
end

info.time = timevector;
info.vars = flux_varnames;
info.ncells = length(fluxnames);
info.nt = length(timevector);
info.flux_out_dir = vic_out_dir;
info.headerlines = headerlines;

info.lon = lon;
info.lat = lat;
% info.lon = flipud(lon);
% info.lat = flipud(lat); % should be in decreasing order
% check_latlon(info.lat, info.lon);

return
