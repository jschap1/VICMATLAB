% Make outputs struct
%
% Makes a structure called OUTPUTS containing key outputs from the VIC
% model run. The structure contains all the information necessary to easily
% make plots of the output variables.

function OUTPUTS = make_outputs_struct(info, wb_avg_ts, wb_avg_map, eb_avg_ts, eb_avg_map)

lon = info.lon;
lat = info.lat;

OUTPUTS = struct();
for p=4:length(info.wbvars)
    OUTPUTS.WB.ts.(info.wbvars{p}) = wb_avg_ts(:,p);
    map1 = xyz2grid(info.lon, info.lat, wb_avg_map(:,p));
    OUTPUTS.WB.maps.(info.wbvars{p}) = map1;   
end

for p=1:length(info.ebvars)
    OUTPUTS.EB.ts.(info.ebvars{p}) = eb_avg_ts(:,p);
    map1 = xyz2grid(lon, lat, eb_avg_map(:,p)); % IRB
    OUTPUTS.EB.maps.(info.ebvars{p}) = map1;
end

OUTPUTS.time = info.time;

OUTPUTS.lat = lat;
OUTPUTS.lon = lon;

% Assign units to output structure
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