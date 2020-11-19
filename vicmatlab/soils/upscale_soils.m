% Upscale soils
% 
% 9/11/2020 JRS
% Upscales the VIC soil parameter file from its current resolution to a
% coarser resolution
%
% Crop and rescale in one step
% INPUTS

function [newsoils] = upscale_soils(soils, newres, oldres)

% Make map
lat_vect = soils(:,3);
lon_vect = soils(:,4);

% Round to the nearest newres increment
newlon1 = floor(min(lon_vect)) + 0.5*newres;
newlon2 = ceil(max(lon_vect)) - 0.5*newres;
newlat1 = floor(min(lat_vect)) + 0.5*newres;
newlat2 = ceil(max(lat_vect)) - 0.5*newres;

newlon = newlon1:newres:newlon2;
newlat = newlat1:newres:newlat2;

% May need to add/subtract half a res here...
R1 = makerefmat(min(lon_vect), min(lat_vect), oldres, oldres); % fine
R2 = makerefmat(newlon(1), newlat(1), newres, newres); % coarse

% Still trying to figure out the best way to upscale... 9/11/2020 JRS

[lon11, lon11] = getlatlon(svar_map, R1);

nvars = size(soils,2);

for k=1:nvars
    svar = soils(:,k);
    svar_map = xyz2grid(lon_vect, lat_vect, svar);
    svar_map2 = svar_map;
    svar_map2(isnan(svar_map2)) = 0;
    % Bringing over a method from an old project
%     addpath('/hdd/SWOTDA/Codes/Met_Forcing_Downscaling/Codes/')
%     [v_interp, finelon, finelat, coarselon, coarselat] = interpolate_merra2(svar_map, R1, R2, 'linear', 0);
    


[Xq,Yq] = ndgrid(lon_vect, lat_vect); % fine
[X,Y] = ndgrid(newlon, newlat); % coarse 

G = griddedInterpolant(X,Y, svar_map, 'linear');
v_interp = G(Xq,Yq);
    
    
    F = griddedInterpolant(svar_map2);
    a = F({newlon,newlat});
end

% Resample map

% Write soil file


return