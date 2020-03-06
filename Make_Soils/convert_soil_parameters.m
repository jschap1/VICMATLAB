% Converts soil parameter file to geotiff files
% 
% INPUTS
% soils = soil parameter file
% setup = argument for get_soil_var_names
% outdir = location to save geotiffs
%
% param = 'ks'
% soilpath = '/Volumes/HD3/VICParametersCONUS';
% soilname = 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt'; 
% infile = fullfile(soilpath, soilname);
% o = 'ks.nc'
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

% tuolumne inputs
% soils = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_L2015/soils_tuolumne_snapped.txt';
% outdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_L2015/soil_plots';
% convert_soil_parameters(soils, '3L-no-org-frost-msds', outdir)

% one-degree tile inputs
% soils = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/soils_snapped.SB';
% outdir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/soil_plots';
% convert_soil_parameters(soils, '3L-no-org-frost-msds', outdir)

% VICGlobal
% soils = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_4/Classic/soils_3L_MERIT.txt';
% outdir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/soil_plots';
% template = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/slope.tif';
% convert_soil_parameters(soils, '3L-no-org-frost-msds', outdir)

% Upper Colorado
% soils = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_vg.txt';
% outdir = '/Volumes/HD4/SWOTDA/Data/Colorado/soils';
% setup = '3L-no-org-frost-msds';
% maskname = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';

% soils = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_soils_L15.txt';
% outdir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/L15/soils';
% mkdir(outdir)
% setup = 'livneh';
% maskname = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_mask.tif';

function [] = convert_soil_parameters(soils, setup, outdir, maskname)

varnames = get_soil_var_names(setup);
[landmask, Rdem] = geotiffread(maskname);
landmask = flipud(landmask);

if ischar(soils)
    soils = load(soils);
    disp('Read soil parameter file')
end

Rdem_mat = georefobj2mat(Rdem, 'LL');
[lon, lat] = pixcenters(Rdem_mat, size(landmask));
landcells = find(landmask == 1);

% lat1 = soils(:,3);
% lon1 = soils(:,4);

% nlat = length(unique(lat));
% nlon = length(unique(lon));

% rasterSize = [nlat, nlon];
% latlim = [min(lat), max(lat)];
% lonlim = [min(lon), max(lon)];

% ulat = unique(lat);
% ulon = unique(lon);

% xres = ulon(2)-ulon(1);
% yres = ulat(2)-ulat(1);

% rasterSize = [nlat, nlon];
% R = georefcells(latlim,lonlim,rasterSize);
% R = georefcells(latlim, lonlim, xres, yres);
    
    for k=5:length(varnames)

        % xyz2grid function from Matlab File Exchange :)
        A = NaN(size(landmask));
%         nodataval = -9999;
%         A = ones(size(landmask))*nodataval;
        
%         A(landcells) = soils(:,k);
        A(landcells) = soils(landcells,k);
% landcells(end) = [];

%         
         
%         A = xyz2grid(lon, lat, soils(:,k)); 
        % 444 rows (latitude) by 922 columns (longitude)
        
        outname = fullfile(outdir, [varnames{k} '.tif']);
        geotiffwrite(outname, flipud(A), Rdem)
        
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