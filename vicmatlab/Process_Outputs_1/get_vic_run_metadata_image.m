function info_image = get_vic_run_metadata_image(fluxfile)

info_image = ncinfo(fluxfile);
info_image.time = ncread(fluxfile, 'time');
info_image.nt = length(info_image.time);

lat1 = ncread(fluxfile, 'lat');
lon1 = ncread(fluxfile, 'lon');

try 
    var_image = ncread(fluxfile, 'OUT_PREC');
end

try 
    var_image = ncread(fluxfile, 'OUT_SWE');
end

[lon_list, lat_list, ~] = grid2xyz(lon1, lat1, var_image(:,:,1));
info_image.lon = lon_list;
info_image.lat = lat_list;

tmp = strsplit(fluxfile, '.');
info_image.out_dir = tmp{1}(1:end-6);

info_image.ncells = length(lon_list);

nc_metadata = ncinfo(fluxfile);

nvars = length(nc_metadata.Variables);
varnames = cell(nvars,1);
for i=1:nvars
    varnames{i} = nc_metadata.Variables(i).Name;
end

info_image.vars = varnames;

return