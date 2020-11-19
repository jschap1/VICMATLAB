% Write soil moisture to soil parameter file
% 
% 8/18/2020 JRS
% Writes initial soil moisture to soil parameter file, keeping 
% the correct spatial arrangement
%
% INPUTS
% spf = name of soil parameter file
% sm = soil moisture data
% lon, lat = column vectors with lon, lat data for each grid cell

function soils = write_soil_moisture_to_spf(spf, sm, lon, lat)

soils = load(spf);

slat = soils(:,3);
slon = soils(:,4);

T = delaunayn([slat, slon]);
ncells = length(slon);
disp(['There are ' num2str(ncells) ' grid cells'])
for j=1:ncells
    k = dsearchn([lat, lon],T,[slat(j) slon(j)]);
%     disp(k)
    soils(j,19) = sm(k,1); % layer 1
    soils(j,20) = sm(k,2); % layer 2
    soils(j,21) = sm(k,3); % layer 3
end

outformat = '3l';
precision = 5;
write_soils(precision, soils, spf, outformat)

return