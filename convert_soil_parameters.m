% Converts soil parameter file to geotiff files
% 
% INPUTS
% infile = soil parameter file
% outfile = savename
% varnames = variable names for the columns of the soil parameter file.
    % Making this an input to keep the function simple

% param = 'ks'
% soilpath = '/Volumes/HD3/VICParametersCONUS';
% soilname = 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt'; 
% infile = fullfile(soilpath, soilname);
% outfile = 'ks.nc'
%
% There is a python tool called tonic that performs a similar function, but
% it is buggy, and I haven't been able to run it on my system 1/16/2018
%
% List of variable names
% allvarnames = ['run_cell', 'gridcel','lat','lon','infilt','Ds','Dsmax', ...
%     'Ws','c','expt','Ksat','phi_s','init_moist','elev','depth', ...
%     'avg_T','dp','bubble','quartz','bulk_density','soil_density', ...
%     'organic', 'bulk_dens_org', 'soil_dens_org', 'off_gmt',...
%     'Wcr_FRACT', 'Wpwp_FRACT', 'rough', 'snow_rough', 'annual_prec', ...
%     'resid_moist', 'fs_active', 'frost_slope', ...
%     'max_snow_distrib_slope', 'July_Tavg'];
%
% Requires mapping toolbox

% varnames = {'run_cell', 'gridcel','lat','lon','infilt','Ds','Dsmax', ...
%     'Ws','c','expt','Ksat','phi_s','init_moist','elev','depth', ...
%     'avg_T','dp','bubble','quartz','bulk_density','soil_density', ...
%     'organic', 'bulk_dens_org', 'soil_dens_org', 'off_gmt',...
%     'Wcr_FRACT', 'Wpwp_FRACT', 'rough', 'snow_rough', 'annual_prec', ...
%     'resid_moist', 'fs_active', 'frost_slope', ...
%     'max_snow_distrib_slope', 'July_Tavg'};

% writecsv_cell('varnames.txt', varnames);

function [] = convert_soil_parameters(infile, varnames)

soils = load(infile);

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
        A = xyz2grid(lon, lat, soils(:,k)); 
        % 444 rows (latitude) by 922 columns (longitude)
        
        outname = [varnames{k} '.tif'];
        geotiffwrite(outname, flipud(A), R)
        disp(['Saved soil parameter data as ', outname])
        
    end

end

% xvals = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(1));
% yvals = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(2));
% imagesc(xvals, yvals, A)

% Scrap ------------------------------------------------------------------

% R.LatitudeLimits(1) % min lat
% R.LatitudeLimits(2) % max lat
% R.LongitudeLimits(1) % min lon
% R.LongitudeLimits(2) % max lon
% 
% [lon, lat, soils(:,k)]