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

% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % no organic, no frozen soil, no July_Tavg, two soil layers
%     'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
%     'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
%     'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
%     'soil_dens1','soil_dens2', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2'};

% writecsv_cell('varnames.txt', varnames);

% infile = '/Volumes/HD3/VICParametersGlobal/vic_params_global_0.5deg/global_soil_param_new';

% infile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/soils_3L_merit_latest.txt';

% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ...
%     'ws','c','expt1','expt2','expt3','ksat1','ksat2','ksat3','phi_s1','phi_s2','phi_s3', ...
%     'init_moist1','init_moist2','init_moist3','elev','depth1','depth2','depth3','avg_T', ...
%     'dp','bubble1','bubble2','bubble3','quartz1','quartz2','quartz3','bulk_dens1','bulk_dens2','bulk_dens3', ...
%     'soil_dens1','soil_dens2','soil_dens3', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wcr_fract3','wpwp_fract1','wpwp_fract2','wpwp_fract3','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2','resid_moist3', 'fs_active'};

% writecsv_cell('varnames.txt', varnames);

% infile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/soils_3L_MERIT_latest.txt';

function [] = convert_soil_parameters(infile, varnames)

soils = load(infile);
disp('loaded soil parameter file')

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