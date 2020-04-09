% Range of soil depths in Livneh et al. VIC setup
% 
% Jan. 10, 2020

[d1, R, lon, lat] = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/depth1.tif');
d2 = flipud(geotiffread('/Volumes/HD3/VICParametersCONUS/soil_data/depth2.tif'));
d3 = flipud(geotiffread('/Volumes/HD3/VICParametersCONUS/soil_data/depth3.tif'));

figure
plotraster(lon, lat, d1, 'Layer 0 thickness (m)','Lon','Lat')

figure
plotraster(lon, lat, d2, 'Layer 1 thickness (m)','Lon','Lat')

figure
plotraster(lon, lat, d3, 'Layer 2 thickness (m)','Lon','Lat')

figure
plotraster(lon, lat, d1+d2+d3, 'Total soil thickness (m)','Lon','Lat')

