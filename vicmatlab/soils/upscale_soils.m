% Upscale soils
% 
% 9/11/2020 JRS
% Upscales the VIC soil parameter file from its current resolution to a
% coarser resolution
%
% INPUTS

function [newsoils, newmask] = upscale_soils(soils, newres, oldres)

% f = newres/oldres; % upscaling (or downscaling) factor
[n, nvars] = size(soils); % n = original number of grid cells

lat_vect = soils(:,3);
lon_vect = soils(:,4);

runcell = xyz2grid(lon_vect, lat_vect, soils(:,1));
figure
plotraster(lon_vect, lat_vect, runcell, 'Runcell');

elev = xyz2grid(lon_vect, lat_vect, soils(:,22));
figure
plotraster(lon_vect, lat_vect, elev, 'Elev');

[runcell2, newlons, newlats] = upscale_raster(runcell, lon_vect, lat_vect, newres, oldres, 'linear');

figure
plotraster(newlons, newlats, runcell2, 'Runcell');

ncells_out = sum(runcell2(:));
disp(['The number of cells in the output is ' num2str(ncells_out)])

R = makerefmat(min(newlons), min(newlats), newres, newres);
extent1 = basin_mask2coordinate_list(runcell2, R);
lat = extent1(:,2); % output coordinates for soil parameter file 
lon = extent1(:,1); 

runcell3 = runcell2(:);
runcell3(runcell3==0) = [];
dat = [lat, lon, runcell3];

runcell4 = xyz2grid(lon, lat, runcell3); % removes empty borders, if they exist
figure,plotraster(lon, lat, runcell4,''); % check

%% Assemble re-scaled soil parameter file

newsoils = zeros(ncells_out, nvars);
newsoils(:,1) = runcell3;
newsoils(:,2) = 1:ncells_out;

for k=5:nvars
    
    invar = xyz2grid(lon_vect, lat_vect, soils(:,k));
    outvar = upscale_raster(invar, lon_vect, lat_vect, newres, oldres, 'linear');
    % outvar = interp2(oldlats, oldlons, invar', newlats, newlons, 'linear')';
    
    % check
%     figure, plotraster(lons, lats, outvar,'');
    outvect = outvar(:);
    
%     lons(runcell2(:)==0) = [];
%     lats(runcell2(:)==0) = [];
    outvect(runcell2(:)==0) = [];
    
%     reconstruct = xyz2grid(lons, lats, outvect);
%     figure
%     plotraster(lons, lats, reconstruct, '');
    
    newsoils(:,k) = outvect;
 
    disp(['Processed soil variable number ' num2str(k) ' of ' num2str(nvars)])
end

lons = zeros(ncells_out,1);
lats = zeros(ncells_out,1);
ind = 1; % brute force fix for lat/lon not matching
for i=1:length(newlons)
    for j=1:length(newlats)
        lons(ind) = newlons(i);
        lats(ind) = newlats(j);
        ind = ind + 1;
    end
end

tmpmask = runcell2(:);
lats(tmpmask==0) = [];
lons(tmpmask==0) = [];
newsoils(:,3) = lats;
newsoils(:,4) = lons;

elev = xyz2grid(lons, lats, newsoils(:,22));
figure
plotraster(lons, lats, elev, 'elev')

newmask = xyz2grid(lons, lats, newsoils(:,1));
figure
plotraster(lons, lats, newmask, 'elev')

return