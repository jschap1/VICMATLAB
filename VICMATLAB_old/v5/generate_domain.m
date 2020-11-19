% Generates a domain file for the VIC 5 image driver
%

%% Create domain file
clc
clear
% fancy but useless
ncdisp(?/Users/dongyueli/Desktop/VIC/src/VIC5.0/VIC_sample_data5.0/image/Stehekin/parameters/domain.stehekin.20151028.nc?);
% simple but on-the-point
ncdisp(?/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc?)
% ##########################################################################
clc
clear
soil=load(?/Users/dongyueli/Dropbox/VICGlobal/data/soils_3L_MERIT.txt?);
[pixel_area,R]=geotiffread(?/Users/dongyueli/Dropbox/VICGlobal/data/pixel_area_km2.tif?);
lat_rec=soil(:,3);      % lat/lon of all pixels, 4141736*1
lon_rec=soil(:,4);
lat=unique(lat_rec);        % unique lat/lon values: 2267*5760
lon=unique(lon_rec);
% get the lon/lat index of each cell in the soil file
for i=1:length(lat_rec)
   lon_ind(i,1)=find(lon==lon_rec(i,1));
   lat_ind(i,1)=find(lat==lat_rec(i,1));
end
% sp_data=dlmread([?disaggregated_forcing/?,file(1).name]);
% time=int32([0:1:(size(sp_data,1)-1)]?);
mask=nan(length(lon),length(lat));
for i=1:length(lat_rec)
   mask(lon_ind(i,1),lat_ind(i,1))=1;
end
frac=nan(length(lon),length(lat));
for i=1:length(lat_rec)
   frac(lon_ind(i,1),lat_ind(i,1))=1;
end
% ATTENTION: THIS PART CAN ONLY BE APPLIED ONCE
pixel_area=flipud(pixel_area);
pixel_area=pixel_area?;
area=nan(length(lon),length(lat));
for i=1:length(lat_rec)
   area(lon_ind(i,1),lat_ind(i,1))=pixel_area(lon_ind(i,1),lat_ind(i,1));
end
area=area*1e6;
% ##########################################################################
nccreate(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?mask?,...
   ?Datatype?,?int32?,...
   ?Dimensions?,{?lon?,length(lon),?lat?,length(lat)},...
         ?Format?,?classic?)
nccreate(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?lon?,...
   ?Datatype?,?double?,...
   ?Dimensions?,{?lon?,length(lon)},...
         ?Format?,?classic?)
nccreate(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?lat?,...
   ?Datatype?,?double?,...
   ?Dimensions?,{?lat?,length(lat)},...
         ?Format?,?classic?)
nccreate(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?frac?,...
   ?Datatype?,?double?,...
   ?Dimensions?,{?lon?,length(lon),?lat?,length(lat)},...
         ?Format?,?classic?)
nccreate(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?area?,...
   ?Datatype?,?double?,...
   ?Dimensions?,{?lon?,length(lon),?lat?,length(lat)},...
         ?Format?,?classic?)
ncwrite(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?mask?,mask);
ncwrite(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?lon?,lon);
ncwrite(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?lat?,lat);
ncwrite(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?frac?,frac);
ncwrite(?/Users/dongyueli/Dropbox/VICGlobal/output/VICGlobal_domain.nc?,?area?,area);