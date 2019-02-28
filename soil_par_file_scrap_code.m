%%
% target_lon = min(lon):0.0625:max(lon);
% target_lat = min(lat):0.0625:max(lat);
% 
% [lons, lats] = ndgrid(min(lon):0.25:max(lon), min(lat):0.25:max(lat));
% [rglons, rglats] = ndgrid(target_lon, target_lat);
% 
% t_sand_rg = myregrid(lons, lats, rglons, rglats, t_sand);

% elev_meta = ncinfo('./Data/JISAO/elev.0.25-deg.nc');
% elev = fliplr(ncread('./Data/JISAO/elev.0.25-deg.nc', 'data')');
% elev = ncread('./Data/JISAO/elev.0.25-deg.nc', 'data');

% lat and lon for elevation grid
% elat = -89.875:0.25:89.875;
% elon = -179.875:0.25:179.875;
% elat = ncread('./Data/JISAO/elev.0.25-deg.nc', 'lat'); 
% elon = ncread('./Data/JISAO/elev.0.25-deg.nc', 'lon');
% elon = 180 - elon;

% Regrid elevation data
[elons, elats] = ndgrid(elon, elat);
elev_rg = myregrid(elons, elats, rglons, rglats, elev);

elev_rg((land_mask==0)) = NaN;

% F = griddedInterpolant(lons, lats, elev');
% elev_rg = F(finelons, finelats);
% tmp = elev_rg(:);
% elev_list = tmp(land_ind);
% clear tmp

% Make sure everything is consistent (it looks good in Panoply, good enough
% for me...)
% figure
% subplot(1,2,1)
% imagesc([lon(1), lon(end)], [lat(1), lat(end)], t_sand)
% subplot(1,2,2)
% imagesc(elev)
% set(gca, 'ydir', 'normal');
% set(gca, 'ydir', 'normal');

%%

%%
addpath('/Users/jschap/Documents/MATLAB/from_file_exchange')

elev = geotiffread();

%% Load elevation data


figure, imagesc(lon, lat, elev), colorbar
set(gca, 'ydir', 'normal');

% Load land mask
% tmp_meta = ncinfo('./Data/JISAO/fractional_land.0.25-deg.nc');
% fract = ncread('./Data/JISAO/fractional_land.0.25-deg.nc', 'data')';
% figure, imagesc(lon, lat, fract)
% set(gca, 'ydir', 'normal');

% Mask the elevation data
% elev(fract==0) = NaN;

elev(~land_ind)
land_ind % from HWSD



%% Regrid to desired resolution

% Use cdo command line tool (if you can get it installed with netcdf4
% support, which I haven't been able to do)
% cdo remapnn,r5760x2880 ../HWSD/HWSD_1247/data/T_SAND.nc4 ../HWSD/HWSD_1247/data/T_SAND_rg.nc4
% Use Matlab griddenInterpolant object

% function rg = myregrid(olons, olats, rglons, rglats, og)
% 
%     F = griddedInterpolant(olons, olats, og);
%     rg = F(rglons, rglats); 
%     
% end

[lons, lats] = ndgrid(min(lon):0.05:max(lon), min(lat):0.05:max(lat));
[rglons, rglats] = ndgrid(min(lon):0.0625:max(lon), min(lat):0.0625:max(lat));
sand_rg = myregrid(lons, lats, rglons, rglats, t_sand);
clay_rg = myregrid(lons, lats, rglons, rglats, t_clay);
land_ind = find(~isnan(sand_rg));


% sand_rg(land_ind);

F = griddedInterpolant(lons, lats, t_sand);

coarsesand = F(finelons, finelats);

figure, imagesc(t_sand)
figure, imagesc(coarsesand)

xyz = raster2xyz(lon, lat, t_clay);
xyz_upscaled = raster2xyz(min(lon):0.0625:max(lon), min(lat):0.0625:max(lat), clay_rg);

% Find the pixels that have data/are not in the ocean or Antarctica
land_ind = find(~isnan(xyz_upscaled(:,3))); 
xyz_upscaled(land_ind,:)
ncells = length(land_ind);

lonlat = xyz(land_ind,1:2);

% lats = meshgrid(lat, lon);
% lons = meshgrid(lon, lat);
% figure, plot3(lons, lats, sand);
% finelat = min(lat):0.0625:max(lat);
% finelon = min(lon):0.0625:max(lon);
% figure, imagesc(finelons, finelats, finesand);

% Regrid elevation data (the current workflow is flawed... doesn't produce
% good soil parameter file)
lat = ncread('/Users/jschap/Downloads/elev.0.25-deg.nc', 'lat');
lon = ncread('/Users/jschap/Downloads/elev.0.25-deg.nc', 'lon');
lon = 180 - lon;
[lons, lats] = ndgrid(min(lon):0.25:max(lon), min(lat):0.25:max(lat));
F = griddedInterpolant(lons, lats, elev');
elev_rg = F(finelons, finelats);
tmp = elev_rg(:);
elev_list = tmp(land_ind);
clear tmp
