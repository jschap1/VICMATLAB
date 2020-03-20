%% Create domain file

clc
clear

% fancy but useless
ncdisp('/Users/dongyueli/Desktop/VIC/src/VIC5.0/VIC_sample_data5.0/image/Stehekin/parameters/domain.stehekin.20151028.nc');
% simple but on-the-point
ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc')

% ##########################################################################

clc
clear

soil=load('/Users/dongyueli/Dropbox/VICGlobal/data/soils_3L_MERIT.txt');
[pixel_area,R]=geotiffread('/Users/dongyueli/Dropbox/VICGlobal/data/pixel_area_km2.tif');

lat_rec=soil(:,3);      % lat/lon of all pixels, 4141736*1    
lon_rec=soil(:,4);

lat=unique(lat_rec);        % unique lat/lon values: 2267*5760
lon=unique(lon_rec);

% get the lon/lat index of each cell in the soil file
for i=1:length(lat_rec)
    lon_ind(i,1)=find(lon==lon_rec(i,1));
    lat_ind(i,1)=find(lat==lat_rec(i,1));
end

% sp_data=dlmread(['disaggregated_forcing/',file(1).name]);
% time=int32([0:1:(size(sp_data,1)-1)]');

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
pixel_area=pixel_area';
area=nan(length(lon),length(lat));
for i=1:length(lat_rec)
    area(lon_ind(i,1),lat_ind(i,1))=pixel_area(lon_ind(i,1),lat_ind(i,1)); 
end
area=area*1e6;
% ##########################################################################


nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','classic')      

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','classic') 
     
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','frac',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','area',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','mask',mask);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','lon',lon);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','lat',lat);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','frac',frac);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_domain.nc','area',area);



%% Create parameter file for CONUS

clc
clear

% ncdisp('/Users/dongyueli/Desktop/VIC/src/VIC5.0/VIC_sample_data5.0/image/Stehekin/parameters/Stehekin_test_params_20160327.nc');
% ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params_conus.nc');

% FillValue_int=-2147483647;
% FillValue_f=9.969209968386869e+36;

% ##########################################################################
soil=dlmread('/Users/dongyueli/Dropbox/VICGlobal/data/soils_3L_MERIT.txt');
snow=dlmread('/Users/dongyueli/Dropbox/VICGlobal/data/snowbands_MERIT.txt');
veg_param=dlmread('/Users/dongyueli/Dropbox/VICGlobal/data/global_vegetation_1_16_IGBP.txt');
% ##########################################################################

% ##########################################################################
% this is to get the parameter file for the entire CONUS
lat=unique(soil(:,3));        % unique lat/lon values: 2267*5760
lon=unique(soil(:,4));

% get the lon/lat index of each cell in the soil file
for i=1:size(soil,1)
    lon_ind(i,1)=find(lon==soil(i,4));
    lat_ind(i,1)=find(lat==soil(i,3));
end
% ##########################################################################

% ##########################################################################
for jj=1:1 % this jj for-loop is only for code folding...

% SOIL FILE PART    
mask=nan(length(lon),length(lat));
for i=1:length(lon_ind)
    mask(lon_ind(i,1),lat_ind(i,1))=1;
end

layer=[1;2;3];

run_cell=nan(length(lon),length(lat));
for i=1:length(lon_ind)
    run_cell(lon_ind(i,1),lat_ind(i,1))=1;
end

gridcell=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    gridcell(lon_ind(i,1),lat_ind(i,1))=soil(i,2);
end

lats=nan(length(lon),length(lat));
for i=1:length(lon)
    lats(i,:)=lat';
end
lats(isnan(mask))=nan;

lons=nan(length(lon),length(lat));
for i=1:length(lat)
    lons(:,i)=lon;
end
lons(isnan(mask))=nan;

infilt=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    infilt(lon_ind(i,1),lat_ind(i,1))=soil(i,5);
end

Ds=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Ds(lon_ind(i,1),lat_ind(i,1))=soil(i,6);
end

Dsmax=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Dsmax(lon_ind(i,1),lat_ind(i,1))=soil(i,7);
end

Ws=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Ws(lon_ind(i,1),lat_ind(i,1))=soil(i,8);
end

c=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    c(lon_ind(i,1),lat_ind(i,1))=soil(i,9);
end

expt=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    expt(lon_ind(i,1),lat_ind(i,1),1)=soil(i,10);
    expt(lon_ind(i,1),lat_ind(i,1),2)=soil(i,11);
    expt(lon_ind(i,1),lat_ind(i,1),3)=soil(i,12);
end

Ksat=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Ksat(lon_ind(i,1),lat_ind(i,1),1)=soil(i,13);
    Ksat(lon_ind(i,1),lat_ind(i,1),2)=soil(i,14);
    Ksat(lon_ind(i,1),lat_ind(i,1),3)=soil(i,15);
end

phi_s=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    phi_s(lon_ind(i,1),lat_ind(i,1),1)=soil(i,16);
    phi_s(lon_ind(i,1),lat_ind(i,1),2)=soil(i,17);
    phi_s(lon_ind(i,1),lat_ind(i,1),3)=soil(i,18);
end

init_moist=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    init_moist(lon_ind(i,1),lat_ind(i,1),1)=soil(i,19);
    init_moist(lon_ind(i,1),lat_ind(i,1),2)=soil(i,20);
    init_moist(lon_ind(i,1),lat_ind(i,1),3)=soil(i,21);
end

elev=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    elev(lon_ind(i,1),lat_ind(i,1))=soil(i,22);
end

depth=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    depth(lon_ind(i,1),lat_ind(i,1),1)=soil(i,23);
    depth(lon_ind(i,1),lat_ind(i,1),2)=soil(i,24);
    depth(lon_ind(i,1),lat_ind(i,1),3)=soil(i,25);
end

avg_T=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    avg_T(lon_ind(i,1),lat_ind(i,1))=soil(i,26);
end

dp=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    dp(lon_ind(i,1),lat_ind(i,1))=soil(i,27);
end

bubble=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    bubble(lon_ind(i,1),lat_ind(i,1),1)=soil(i,28);
    bubble(lon_ind(i,1),lat_ind(i,1),2)=soil(i,29);
    bubble(lon_ind(i,1),lat_ind(i,1),3)=soil(i,30);
end

quartz=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    quartz(lon_ind(i,1),lat_ind(i,1),1)=soil(i,31);
    quartz(lon_ind(i,1),lat_ind(i,1),2)=soil(i,32);
    quartz(lon_ind(i,1),lat_ind(i,1),3)=soil(i,33);
end

bulk_density=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    bulk_density(lon_ind(i,1),lat_ind(i,1),1)=soil(i,34);
    bulk_density(lon_ind(i,1),lat_ind(i,1),2)=soil(i,35);
    bulk_density(lon_ind(i,1),lat_ind(i,1),3)=soil(i,36);
end

soil_density=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    soil_density(lon_ind(i,1),lat_ind(i,1),1)=soil(i,37);
    soil_density(lon_ind(i,1),lat_ind(i,1),2)=soil(i,38);
    soil_density(lon_ind(i,1),lat_ind(i,1),3)=soil(i,39);
end

off_gmt=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    off_gmt(lon_ind(i,1),lat_ind(i,1))=soil(i,40);
end

Wcr_FRACT=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),1)=soil(i,41);
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),2)=soil(i,42);
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),3)=soil(i,43);
end

Wpwp_FRACT=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),1)=soil(i,44);
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),2)=soil(i,45);
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),3)=soil(i,46);
end

rough=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    rough(lon_ind(i,1),lat_ind(i,1))=soil(i,47);
end

snow_rough=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    snow_rough(lon_ind(i,1),lat_ind(i,1))=soil(i,48);
end

annual_prec=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    annual_prec(lon_ind(i,1),lat_ind(i,1))=soil(i,49);
end

resid_moist=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    resid_moist(lon_ind(i,1),lat_ind(i,1),1)=soil(i,50);
    resid_moist(lon_ind(i,1),lat_ind(i,1),2)=soil(i,51);
    resid_moist(lon_ind(i,1),lat_ind(i,1),3)=soil(i,52);
end

fs_active=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    fs_active(lon_ind(i,1),lat_ind(i,1))=soil(i,53);
end

July_Tavg=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    July_Tavg(lon_ind(i,1),lat_ind(i,1))=soil(i,54);
end

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','layer',...
    'Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},...
          'Format','netcdf4_classic')
          
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','run_cell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','gridcell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lats',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lons',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','infilt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ds',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Dsmax',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ws',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','c',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','expt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ksat',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','phi_s',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','init_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','elev',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','avg_T',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','dp',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')      

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','bubble',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','quartz',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','bulk_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','soil_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','off_gmt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Wcr_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Wpwp_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','snow_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')      

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','annual_prec',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','resid_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','fs_active',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','July_Tavg',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lat',lat);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lon',lon);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','mask',mask);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','layer',layer);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','run_cell',run_cell);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','gridcell',gridcell);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lats',lats);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','lons',lons);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','infilt',infilt);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ds',Ds);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Dsmax',Dsmax);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ws',Ws);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','c',c);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','expt',expt);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Ksat',Ksat);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','phi_s',phi_s);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','init_moist',init_moist);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','elev',elev);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','depth',depth);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','avg_T',avg_T);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','dp',dp);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','bubble',bubble);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','quartz',quartz);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','bulk_density',bulk_density);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','soil_density',soil_density);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','off_gmt',off_gmt);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Wcr_FRACT',Wcr_FRACT);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Wpwp_FRACT',Wpwp_FRACT);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rough',rough);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','snow_rough',snow_rough);   
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','annual_prec',annual_prec);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','resid_moist',resid_moist);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','fs_active',fs_active);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','July_Tavg',July_Tavg);

clear July_Tavg fs_active resid_moist annual_prec snow_rough rough Wpwp_FRACT 
clear soil_density bulk_density quartz bubble dp avg_T depth elev init_moist gridcell
clear mask run_cell Wcr_FRACT off_gmt Dsmax Ds lons phi_s Ksat expt c Ws infilt lats

end


for jj=1:1 % this jj for-loop is only for code folding...

% SNOW FILE PART
snow_band=[1;2;3;4;5];

cellnum=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    cellnum(lon_ind(i,1),lat_ind(i,1))=snow(i,1);
end

AreaFract=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    AreaFract(lon_ind(i,1),lat_ind(i,1),1)=snow(i,2);
    AreaFract(lon_ind(i,1),lat_ind(i,1),2)=snow(i,3);
    AreaFract(lon_ind(i,1),lat_ind(i,1),3)=snow(i,4);
    AreaFract(lon_ind(i,1),lat_ind(i,1),4)=snow(i,5);
    AreaFract(lon_ind(i,1),lat_ind(i,1),5)=snow(i,6);
end

elevation=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    elevation(lon_ind(i,1),lat_ind(i,1),1)=snow(i,7);
    elevation(lon_ind(i,1),lat_ind(i,1),2)=snow(i,8);
    elevation(lon_ind(i,1),lat_ind(i,1),3)=snow(i,9);
    elevation(lon_ind(i,1),lat_ind(i,1),4)=snow(i,10);
    elevation(lon_ind(i,1),lat_ind(i,1),5)=snow(i,11);
end

Pfactor=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    Pfactor(lon_ind(i,1),lat_ind(i,1),1)=snow(i,12);
    Pfactor(lon_ind(i,1),lat_ind(i,1),2)=snow(i,13);
    Pfactor(lon_ind(i,1),lat_ind(i,1),3)=snow(i,14);
    Pfactor(lon_ind(i,1),lat_ind(i,1),4)=snow(i,15);
    Pfactor(lon_ind(i,1),lat_ind(i,1),5)=snow(i,16);
end

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','snow_band',...
    'Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','cellnum',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','AreaFract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','elevation',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','netcdf4_classic')
     
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Pfactor',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','netcdf4_classic')

ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','snow_band',snow_band);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','cellnum',cellnum);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','AreaFract',AreaFract);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','elevation',elevation);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Pfactor',Pfactor);

clear cellnum AreaFract elevation Pfactor

end


for jj=1:1 % this jj for-loop is only for code folding...



% VEG FILE AND VEG LIB PART
nh_veg=load('/Users/dongyueli/Dropbox/VICGlobal/data/veglib_nh_nohead.txt');
sh_veg=load('/Users/dongyueli/Dropbox/VICGlobal/data/veglib_sh_nohead.txt');
nh_veg_lib=[nh_veg(:,1:16) nh_veg(:,29:69)];
sh_veg_lib=[sh_veg(:,1:16) sh_veg(:,29:69)];
nh_LAI=nh_veg(:,5:16);
sh_LAI=sh_veg(:,5:16);
nh_fcan=nh_veg(:,17:28);
sh_fcan=sh_veg(:,17:28);

veg_class=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17];

veg_descr=char({'Open Water','Evergreen Needleleaf','Evergreen Broadleaf','Deciduous Needleleaf',...
    'Deciduous Broadleaf','Mixed Forest','Closed Shrub','Open Shrub','Savannas','Woody Savannas',...
    'Grasslands','Permanent Wetlands','Crop land','Urban','Cropland NAtural Veg Mosaic',...
    'SnowIce','Barren'});

root_zone=[1;2;3];

month=[1;2;3;4;5;6;7;8;9;10;11;12];

line_no=find(veg_param(:,3)==0);
veg_param_line=veg_param(line_no,:);

Nveg=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Nveg(lon_ind(i,1),lat_ind(i,1))=veg_param_line(i,2);
end

Cv=nan(length(lon),length(lat),length(veg_class));
root_depth=nan(length(lon),length(lat),length(root_zone),length(veg_class));
root_fract=nan(length(lon),length(lat),length(root_zone),length(veg_class));
LAI=nan(length(lon),length(lat),length(month),length(veg_class));

for i=1:length(lon_ind)
    tile_no=veg_param_line(i,2);
    pix_start_line=line_no(i,1);
    for j=1:tile_no
        Cv(lon_ind(i,1),lat_ind(i,1),veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,2);
        root_depth(lon_ind(i,1),lat_ind(i,1),1,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,3);
        root_fract(lon_ind(i,1),lat_ind(i,1),1,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,4);
        root_depth(lon_ind(i,1),lat_ind(i,1),2,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,5);
        root_fract(lon_ind(i,1),lat_ind(i,1),2,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,6);
        root_depth(lon_ind(i,1),lat_ind(i,1),3,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,7);
        root_fract(lon_ind(i,1),lat_ind(i,1),3,veg_param(pix_start_line+1*(j-1)+1,1))...
            =veg_param(pix_start_line+1*(j-1)+1,8);        
        if soil(i,3)<=0
            LAI(lon_ind(i,1),lat_ind(i,1),1:12,veg_param(pix_start_line+1*(j-1)+1,1))...
                =sh_LAI(veg_param(pix_start_line+1*(j-1)+1,1),1:12);        
        else
            LAI(lon_ind(i,1),lat_ind(i,1),1:12,veg_param(pix_start_line+1*(j-1)+1,1))...
                =nh_LAI(veg_param(pix_start_line+1*(j-1)+1,1),1:12);
        end
    end   
end

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_class',...
    'Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_descr',...
    'Datatype','char',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_zone',...
    'Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},...
          'Format','netcdf4_classic')
      
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','month',...
    'Datatype','int32',...
    'Dimensions',{'month',length(month)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Nveg',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Cv',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_fract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','LAI',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')

ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_class',veg_class);
disp('1')
% ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_descr',veg_descr);
disp('2')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_zone',root_zone);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','month',month);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Nveg',Nveg);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','Cv',Cv);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_depth',root_depth);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','root_fract',root_fract);
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','LAI',LAI);

clear LAI root_depth root_fract


overstory=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        overstory(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,2);
    else
        overstory(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,2);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','overstory',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','overstory',overstory);
clear overstory


rarc=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        rarc(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,3);
    else
        rarc(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,3);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rarc',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rarc',rarc);
clear rarc

rmin=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)    
    if soil(i,3)<=0
        rmin(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,4);
    else
        rmin(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,4);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rmin',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rmin',rmin);
clear rmin

wind_h=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)   
    if soil(i,3)<=0
        wind_h(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,53);
    else
        wind_h(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,53); 
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','wind_h',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','wind_h',wind_h);
clear wind_h

RGL=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        RGL(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,54);
    else
        RGL(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,54);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','RGL',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','RGL',RGL);
clear RGL

rad_atten=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        rad_atten(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,55);
    else
        rad_atten(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,55); 
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rad_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','rad_atten',rad_atten);    
clear rad_atten

wind_atten=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)    
    if soil(i,3)<=0
        wind_atten(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,56);
    else
        wind_atten(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,56);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','wind_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','wind_atten',wind_atten);
clear wind_atten
      
trunk_ratio=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)  
    if soil(i,3)<=0
        trunk_ratio(lon_ind(i,1),lat_ind(i,1),:)=sh_veg_lib(:,57);
    else
        trunk_ratio(lon_ind(i,1),lat_ind(i,1),:)=nh_veg_lib(:,57);
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','trunk_ratio',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','trunk_ratio',trunk_ratio);
clear trunk_ratio      

albedo=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        for j=1:length(veg_class)
            albedo(lon_ind(i,1),lat_ind(i,1),:,j)=sh_veg_lib(j,17:28);
        end
    else
        for j=1:length(veg_class)
            albedo(lon_ind(i,1),lat_ind(i,1),:,j)=nh_veg_lib(j,17:28);
        end 
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','albedo',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','albedo',albedo);
clear albedo

veg_rough=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        for j=1:length(veg_class)
            veg_rough(lon_ind(i,1),lat_ind(i,1),:,j)=sh_veg_lib(j,29:40);
        end
    else
        for j=1:length(veg_class)
            veg_rough(lon_ind(i,1),lat_ind(i,1),:,j)=nh_veg_lib(j,29:40);
        end
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','veg_rough',veg_rough);
clear veg_rough      

displacement=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        for j=1:length(veg_class)
            displacement(lon_ind(i,1),lat_ind(i,1),:,j)=sh_veg_lib(j,41:52);
        end
    else
        for j=1:length(veg_class)
            displacement(lon_ind(i,1),lat_ind(i,1),:,j)=nh_veg_lib(j,41:52);
        end                
    end
end
nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','displacement',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')  
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','displacement',displacement);
clear displacement      

fcanopy=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if soil(i,3)<=0
        for j=1:length(veg_class)
            fcanopy(lon_ind(i,1),lat_ind(i,1),:,j)=sh_fcan(j,1:12);
        end
    else
        for j=1:length(veg_class)
            fcanopy(lon_ind(i,1),lat_ind(i,1),:,j)=nh_fcan(j,1:12);
        end 
    end
end

nccreate('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','fcanopy',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','netcdf4_classic')  
ncwrite('/Volumes/Disk4/VICGlobal/output/VICGlobal_params.nc','fcanopy',fcanopy);
clear fcanopy

end




