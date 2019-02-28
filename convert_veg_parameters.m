% Converts vegetation parameter file to a geotiff file
% 
% INPUTS
% vegfile = vegetation parameter file
% soilfile = soil parameter file
% outfile = savename
% processed vegetation parameter file from load_veg_parameters.m
%
% Requires mapping toolbox
% soilfile = 'global_soil_param_new';
% vegfile = 'global_veg_param_new';



soilfile = 'global_soil_param_new';
soils = load(soilfile);

lat = soils(:,3);
lon = soils(:,4);

nlat = length(unique(lat));
nlon = length(unique(lon));

rasterSize = [nlat, nlon];
latlim = [min(lat), max(lat)];
lonlim = [min(lon), max(lon)];

R = georefcells(latlim,lonlim,rasterSize);
    
    for k=5:length(varnames)

        % xyz2grid function from Matlab File Exchange :)
%         A = xyz2grid(lon, lat, soils(:,k)); 
        A = xyz2grid(lon, lat, nveg'); 
        % 444 rows (latitude) by 922 columns (longitude)
        
        outname = [varnames{k} '.tif'];
        geotiffwrite(outname, flipud(A), R)
        disp(['Saved soil parameter data as ', outname])
        
    end


% xvals = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(1));
% yvals = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(2));
% imagesc(xvals, yvals, A)