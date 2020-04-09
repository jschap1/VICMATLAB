% Snap soils
%
% It is necessary to have matching coordinates between the forcing files
% and the soil parameter file for the VIC model to run. This function 
% finds the closest MERRA-2 forcing file for each row in the soil parameter file
% and changes the coordinates of the soil parameter file to match. 
%
% INPUTS
% soils = soil parameter file
% latlon = coordinates from the MERRA-2 grid/forcing files
% precision = number of decimal places for coordinates in soil parameter
% file
% outname = filename to save the soil parameter file
% outformat = output format for the soil parameter file
%
% OUTPUTS
% newcoords = lat and lon of the new coordinates for the soil parameter
% file
% oldcoords = original coordinates of the soil parameter file
%
% TODO
% Integrate this into clip_soils, get rid of the unnecessary name changing

function [soils, newcoords, oldcoords] = snap_soils(soils, latlon, precision, outname, outformat)

slat = soils(:,3);
slon = soils(:,4);
slatlon = [slat, slon];

flat = latlon(:,1);
flon = latlon(:,2);
flatlon = latlon;

T = delaunayn(flatlon);
k = dsearchn(flatlon,T,slatlon);

soils(:,3) = flat(k);
soils(:,4) = flon(k);

oldcoords = [slat, slon];
newcoords = latlon;

write_soils(precision, soils, outname, outformat)
disp(['Snapped soil parameter file written to ' outname])

return