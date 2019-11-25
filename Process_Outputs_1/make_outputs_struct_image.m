function OUTPUTS = make_outputs_struct_image(info, avg_ts, avg_map)

OUTPUTS.avg_ts = avg_ts;
OUTPUTS.avg_map = avg_map;
OUTPUTS.time = info.time;
OUTPUTS.lon = info.lon;
OUTPUTS.lat = info.lat;

% TODO add units
OUTPUTS.units = {''};

return