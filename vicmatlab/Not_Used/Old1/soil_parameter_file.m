% This is to make the soil parameter file for the Kaweah basin VIC model

% Each grid cell in the model domain gets its own row in the soil
% parameter file. There are 103 1/16 degree grid cells in the Kaweah basin,
% so there should be 103 rows in this soil parameter file.
%
% See list of outputs, at bottom.
%
% For more info on the soil parameter file, see:
% http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/PrepSoilParam.shtml
% http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/SoilParam.shtml









% Attempting to read in the sample soil parameter file from the VIC
% website...
% fname = '/Users/jschapMac/Documents/HydrologyData/soil.param.sample.txt';
% fID = fopen(fname);
% sample = textscan(fID, formatspec, );
























% ------------------------------------------------------------------------

% OUTPUTS
% run_cell = 1 if running cell, 0 otherwise
% gridcel
% lat
% lon
% infilt
% Ds
% Dsmax
% Ws
% c
% expt
% Ksat
% phi_s
% init_moist
% elev
% depth
% avg_T
% dp
% bubble
% quartz
% bulk_density
% soil_density
% organic (optional)
% bulk_dens_org (optional)
% soil_dens_org (optional)
% off_gmt
% Wcr_FRACT
% Wpwp_FRACT
% rough
% snow_rough
% annual_prec
% resid_moist
% fs_active
% frost_slope (optional)
% max_snow_distrib_slope (optional)
% July_Tavg (optional)