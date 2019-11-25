% Make outputs struct
%
% Makes a structure called OUTPUTS containing key outputs from the VIC
% model run. The structure contains all the information necessary to easily
% make plots of the output variables.

function OUTPUTS = make_outputs_struct(info, wb_avg_ts, wb_avg_map, eb_avg_ts, eb_avg_map)

% time series

% [lon, lat] = get_coordinates_from_VIC_file(info.eb_out_dir, 'eb_', 'fluxes');
[lon, lat] = get_coordinates_from_VIC_file(info.eb_out_dir, 'eb_', 'forcings');

OUTPUTS = struct();
for p=4:length(info.wbvars)
    OUTPUTS.WB.ts.(info.wbvars{p}) = wb_avg_ts(:,p);
%     map1 = fliplr(xyz2grid(info.lon, info.lat, wb_avg_map(:,p)));
    map1 = xyz2grid(lon, lat, wb_avg_map(:,p)); % IRB
%     map1 = fliplr(xyz2grid(lon, lat, wb_avg_map(:,p))); % UMRB
    OUTPUTS.WB.maps.(info.wbvars{p}) = map1;  
end

% A = load('/Volumes/HD4/SWOTDA/Data/UMRB/Classic_Livneh_met_L15/Raw/Processed/readVIC_outputs.mat');
% A.ds.Files(1,:)
% 
% info.lon(1)
% info.lat(1)
% 
% figure, plotraster(lon, lat, map1, 'Evaporation (mm)', '','')

for p=1:length(info.ebvars)
    OUTPUTS.EB.ts.(info.ebvars{p}) = eb_avg_ts(:,p);
%     map1 = fliplr(xyz2grid(info.lon, info.lat, eb_avg_map(:,p)));
    map1 = xyz2grid(lon, lat, eb_avg_map(:,p)); % IRB
%     map1 = fliplr(xyz2grid(lon, lat, eb_avg_map(:,p))); % UMRB
    OUTPUTS.EB.maps.(info.ebvars{p}) = map1;
end

OUTPUTS.time = info.time;
% OUTPUTS.pixlat = lat;
% OUTPUTS.pixlon = lon;
% OUTPUTS.lat = info.lat;
% OUTPUTS.lon = info.lon;

OUTPUTS.lat = lat;
OUTPUTS.lon = lon;

% Need to check these
OUTPUTS.units.EVAP = 'mm';
OUTPUTS.units.RUNOFF = 'mm';
OUTPUTS.units.BASEFLOW = 'mm';
OUTPUTS.units.EVAP_CANOP = 'mm';
OUTPUTS.units.TRANSP_VEG = 'mm';
OUTPUTS.units.EVAP_BARE = 'mm';
OUTPUTS.units.SUB_CANOP = 'mm';
OUTPUTS.units.SUB_SNOW = 'mm';
OUTPUTS.units.PET = 'mm';
OUTPUTS.units.PREC = 'mm';
OUTPUTS.units.RAINF = 'mm';
OUTPUTS.units.SNOWF = 'mm';
OUTPUTS.units.WATER_ERROR = 'mm';
OUTPUTS.units.WDEW = 'mm';
OUTPUTS.units.SOIL_LIQ = 'mm';
OUTPUTS.units.SOIL_MOIST = 'mm';
OUTPUTS.units.SNOW_COVER = 'mm';
OUTPUTS.units.SNOW_DEPTH = 'mm';
OUTPUTS.units.SNOW_CANOPY = 'mm';
OUTPUTS.units.SWE = 'mm';

return

%% Comments

% map1 = xyz2grid(info.lon, info.lat, eb_avg_map(:,15));
% figure, plotraster(info.lon, info.lat, fliplr(map1), '','','')
%
% lon decreasing, lat increasing
% 
% map1 = xyz2grid(flipud(lon), lat, eb_avg_map(:,15));
% figure, plotraster(flipud(lon), lat, map1, '','','')
% % lon increasing, lat increasing
% 
% map1 = xyz2grid(flipud(lon), flipud(lat), eb_avg_map(:,15));
% figure, plotraster(flipud(lon), flipud(lat), map1, '','','')
% % lon increasing, lat decreasing
%
% not sure if necessary to flip around... might be a problem here...
% figure, plotraster(lon, lat, fliplr(flipud(map1)), '','','')
%
% Yeah, comparing with the Tuolumne NetCDF plots from the image driver,
% these maps are definitely rotated the wrong way. 
% 1. Find out if this is a general or specific problem (does it apply for
% the one degree test tile?)
% 2. Find out where the problem is occurring and fix it
% For the one-degree test, I did not encounter this problem.
% The soil parameter plots are upside down.
% This implies there is something wrong with the global soil parameter file
% because the subsetting works fine.
% No, actually the L2015 plots are upside down, too. 
%
% OK, so the solution is that when I read in geotiffs, I need to use
% flipud, always.

% figure
% plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.EB.maps.OUT_AIR_TEMP, 'Average air temperature (deg C)', 'Lon', 'Lat')
% 
% figure
% jsplot(OUTPUTS.time, OUTPUTS.EB.ts.OUT_AIR_TEMP, 'VAR', 'Time', 'VAR (units)', 14)

% jsplot(OUTPUTS.time, OUTPUTS.WB.ts.OUT_PREC, 'Precipitation', 'Time', 'Precip. (mm/hr)', 14)
% plotraster(OUTPUTS.lon, OUTPUTS.lat, OUTPUTS.WB.maps.OUT_PREC, 'Average Precipitation (mm/hr)', 'Lon', 'Lat')

% wb_avg_map_table = array2table(wb_avg_map);
% wb_avg_map_table.Properties.VariableNames = info.wbvars;
% % make maps for each variable using xyz2grid
% precip = xyz2grid(info.lon, info.lat, info.eb.avgmap.OUT_PREC);
% 
% wb_avg_ts_table = array2table(wb_avg_ts);
% wb_avg_ts_table.Properties.VariableNames = info.wbvars;

% load('/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/Processed/vic_outputs_summarized_1980-2018_daily_WB.mat')

% load /Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_VICGlobal/Processed_EB_SB/vic_outputs_summarized_1999_daily.mat

% results_dir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/Processed_WB';
% load(fullfile(results_dir, 'vic_outputs_summarized_1999-2001_daily.mat'), 'OUTPUTS');
