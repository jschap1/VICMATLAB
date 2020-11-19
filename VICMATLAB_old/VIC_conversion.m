%% Create domain file

clc
clear

% ##########################################################################
file=dir('forcing_474/data*'); % this is the forcing data directory
for i=1:length(file)
    lat_rec(i,1)=str2num(file(i).name(6:13));    
    lon_rec(i,1)=str2num(file(i).name(15:24));   
end

min_lat=min(lat_rec);
max_lat=max(lat_rec);
min_lon=min(lon_rec);
max_lon=max(lon_rec);

lat=([min_lat:0.0625:max_lat])';
lon=([min_lon:0.0625:max_lon])';

% get the lon/lat index of each cell in the soil file
for i=1:length(file)
    lon_ind(i,1)=find(lon==lon_rec(i,1));
    lat_ind(i,1)=find(lat==lat_rec(i,1));
end

% sp_data=dlmread(['disaggregated_forcing/',file(1).name]);
% time=int32([0:1:(size(sp_data,1)-1)]');

mask=nan(length(lon),length(lat));
for i=1:length(file)
    mask(lon_ind(i,1),lat_ind(i,1))=1; 
end

frac=nan(length(lon),length(lat));
for i=1:length(file)
    frac(lon_ind(i,1),lat_ind(i,1))=1; 
end

area=nan(length(lon),length(lat));
for i=1:length(file)
    area(lon_ind(i,1),lat_ind(i,1))=37300000.0; 
end
% ##########################################################################


nccreate('Shasta_VIC5_domain.nc','mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')


nccreate('Shasta_VIC5_domain.nc','lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','classic')      

nccreate('Shasta_VIC5_domain.nc','lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','classic') 
     
nccreate('Shasta_VIC5_domain.nc','frac',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate('Shasta_VIC5_domain.nc','area',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')


ncwrite('Shasta_VIC5_domain.nc','mask',mask);
ncwrite('Shasta_VIC5_domain.nc','lon',lon);
ncwrite('Shasta_VIC5_domain.nc','lat',lat);
ncwrite('Shasta_VIC5_domain.nc','frac',frac);
ncwrite('Shasta_VIC5_domain.nc','area',area);





ncdisp('./shasta_vic5/Shasta_VIC5_domain.nc')

a=[];
a=ncread('./shasta_vic5/Shasta_VIC5_domain.nc','lon');




% working Sierra Nevada version
ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc')

a=[];
a=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc','mask');




%% Create parameter file for the domain (same extent with domain file)

clc
clear

% FillValue_int=-2147483647;
% FillValue_f=9.969209968386869e+36;

% ##########################################################################
load('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/vegetation_library.mat');
soil=dlmread('/Users/dongyueli/Desktop/VIC/setup/Kostas_sierra_vic/input/vic.soil.0625.new.cal.adj.conus.plus.crb.can');
soil=soil(1:208326,:); % from 208327 to the end of the file is in Canada
snow=dlmread('/Users/dongyueli/Desktop/VIC/setup/Kostas_sierra_vic/input/vic.snow.0625.new.cal.adj.can.5bands');
snow=snow(1:208326,:);
veg_param=dlmread('/Users/dongyueli/Desktop/VIC/setup/Kostas_sierra_vic/input/vic.veg.0625.new.cal.adj.can');
% ##########################################################################

% ##########################################################################
file=dir('forcing_474/data*'); % this is the forcing data directory
for i=1:length(file)
    lat_rec(i,1)=str2num(file(i).name(6:13));    
    lon_rec(i,1)=str2num(file(i).name(15:24));   
end

min_lat=min(lat_rec); % domain boundary
max_lat=max(lat_rec);
min_lon=min(lon_rec);
max_lon=max(lon_rec);
lat=([min_lat:0.0625:max_lat])'; % domain lat/lon
lon=([min_lon:0.0625:max_lon])';
clear min* max*

for i=1:length(file)
    lon_ind(i,1)=find(lon==lon_rec(i,1)); % position of running pix in DOMAIN
    lat_ind(i,1)=find(lat==lat_rec(i,1));
end

for i=1:length(file)
    k=find(soil(:,4)==lon_rec(i,1)); % position of running pix in SOIL
    t=find(soil(k,3)==lat_rec(i,1));
    soil_ind(i,1)=k(t);
end
clear k t
% ##########################################################################

% ##########################################################################
for jj=1:1 % this jj for-loop is only for code folding...

% SOIL FILE
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
    gridcell(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),2);
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
    infilt(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),5);
end

Ds=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Ds(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),6);
end

Dsmax=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Dsmax(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),7);
end

Ws=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    Ws(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),8);
end

c=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    c(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),9);
end

expt=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    expt(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),10);
    expt(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),11);
    expt(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),12);
end

Ksat=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Ksat(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),13);
    Ksat(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),14);
    Ksat(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),15);
end

phi_s=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    phi_s(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),16);
    phi_s(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),17);
    phi_s(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),18);
end

init_moist=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    init_moist(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),19);
    init_moist(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),20);
    init_moist(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),21);
end

elev=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    elev(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),22);
end

depth=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    depth(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),23);
    depth(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),24);
    depth(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),25);
end

avg_T=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    avg_T(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),26);
end

dp=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    dp(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),27);
end

bubble=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    bubble(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),28);
    bubble(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),29);
    bubble(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),30);
end

quartz=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    quartz(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),31);
    quartz(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),32);
    quartz(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),33);
end

bulk_density=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    bulk_density(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),34);
    bulk_density(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),35);
    bulk_density(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),36);
end

soil_density=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    soil_density(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),37);
    soil_density(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),38);
    soil_density(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),39);
end

off_gmt=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    off_gmt(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),40);
end

Wcr_FRACT=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),41);
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),42);
    Wcr_FRACT(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),43);
end

Wpwp_FRACT=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),44);
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),45);
    Wpwp_FRACT(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),46);
end

rough=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    rough(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),47);
end

snow_rough=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    snow_rough(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),48);
end

annual_prec=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    annual_prec(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),49);
end

resid_moist=nan(length(lon),length(lat),length(layer));
for i=1:length(lon_ind)    
    resid_moist(lon_ind(i,1),lat_ind(i,1),1)=soil(soil_ind(i),50);
    resid_moist(lon_ind(i,1),lat_ind(i,1),2)=soil(soil_ind(i),51);
    resid_moist(lon_ind(i,1),lat_ind(i,1),3)=soil(soil_ind(i),52);
end

fs_active=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    fs_active(lon_ind(i,1),lat_ind(i,1))=soil(soil_ind(i),53);
end


% SNOW FILE
snow_band=[1;2;3;4;5];

cellnum=nan(length(lon),length(lat));
for i=1:length(lon_ind)    
    cellnum(lon_ind(i,1),lat_ind(i,1))=snow(soil_ind(i),1);
end

AreaFract=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    AreaFract(lon_ind(i,1),lat_ind(i,1),1)=snow(soil_ind(i),2);
    AreaFract(lon_ind(i,1),lat_ind(i,1),2)=snow(soil_ind(i),3);
    AreaFract(lon_ind(i,1),lat_ind(i,1),3)=snow(soil_ind(i),4);
    AreaFract(lon_ind(i,1),lat_ind(i,1),4)=snow(soil_ind(i),5);
    AreaFract(lon_ind(i,1),lat_ind(i,1),5)=snow(soil_ind(i),6);
end

elevation=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    elevation(lon_ind(i,1),lat_ind(i,1),1)=snow(soil_ind(i),7);
    elevation(lon_ind(i,1),lat_ind(i,1),2)=snow(soil_ind(i),8);
    elevation(lon_ind(i,1),lat_ind(i,1),3)=snow(soil_ind(i),9);
    elevation(lon_ind(i,1),lat_ind(i,1),4)=snow(soil_ind(i),10);
    elevation(lon_ind(i,1),lat_ind(i,1),5)=snow(soil_ind(i),11);
end

Pfactor=nan(length(lon),length(lat),length(snow_band));
for i=1:length(lon_ind)    
    Pfactor(lon_ind(i,1),lat_ind(i,1),1)=snow(soil_ind(i),12);
    Pfactor(lon_ind(i,1),lat_ind(i,1),2)=snow(soil_ind(i),13);
    Pfactor(lon_ind(i,1),lat_ind(i,1),3)=snow(soil_ind(i),14);
    Pfactor(lon_ind(i,1),lat_ind(i,1),4)=snow(soil_ind(i),15);
    Pfactor(lon_ind(i,1),lat_ind(i,1),5)=snow(soil_ind(i),16);
end

% VEG FILE
veg_class=[1;2;3;4;5;6;7;8;9;10;11;12];

veg_descr=char({'Evergreen Needleleaf','Evergreen Broadleaf','Deciduous Needleleaf',...
    'Deciduous Broadleaf','Mixed Cover','Woodland','Wooded Grassland',...
    'Closed Shrubland','Open Shrubland','Grasslands','Crop land','Bare Soil'});

root_zone=[1;2;3];

month=[1;2;3;4;5;6;7;8;9;10;11;12];

% grab the pixel lines in veg_param file
line_no=find(veg_param(:,3)==0);
veg_param_line=veg_param(line_no,1:2);

Nveg=nan(length(lon),length(lat));
for i=1:length(lon_ind)
    Nveg(lon_ind(i,1),lat_ind(i,1))=veg_param_line(soil_ind(i),2);
end

Cv=nan(length(lon),length(lat),length(veg_class));
root_depth=nan(length(lon),length(lat),length(root_zone),length(veg_class));
root_fract=nan(length(lon),length(lat),length(root_zone),length(veg_class));
LAI=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    
    tile_no=veg_param_line(soil_ind(i),2);
    pix_start_line=line_no(soil_ind(i),1);
    
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        Cv(lon_ind(i,1),lat_ind(i,1),:)=0;
        root_fract(lon_ind(i,1),lat_ind(i,1),:,:)=0;
        root_depth(lon_ind(i,1),lat_ind(i,1),:,:)=0;
        LAI(lon_ind(i,1),lat_ind(i,1),:,:)=0;
    end
    
    for j=1:tile_no
        Cv(lon_ind(i,1),lat_ind(i,1),veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,2);
        root_depth(lon_ind(i,1),lat_ind(i,1),1,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,3);
        root_fract(lon_ind(i,1),lat_ind(i,1),1,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,4);
        root_depth(lon_ind(i,1),lat_ind(i,1),2,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,5);
        root_fract(lon_ind(i,1),lat_ind(i,1),2,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,6);
        root_depth(lon_ind(i,1),lat_ind(i,1),3,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,7);
        root_fract(lon_ind(i,1),lat_ind(i,1),3,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+1,8);
        
        LAI(lon_ind(i,1),lat_ind(i,1),1:12,veg_param(pix_start_line+2*(j-1)+1,1))...
            =veg_param(pix_start_line+2*(j-1)+2,1:12);        
    end
    
    if nansum(Cv(lon_ind(i,1),lat_ind(i,1),1:11))<1
        Cv(lon_ind(i,1),lat_ind(i,1),12)=1-nansum(Cv(lon_ind(i,1),lat_ind(i,1),1:11));
    end
    
end

overstory=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)    
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        overstory(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    overstory(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,2);
end

rarc=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        rarc(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    rarc(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,3);
    rarc(lon_ind(i,1),lat_ind(i,1),12)=100;
end

rmin=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        rmin(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    rmin(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,4);
end

wind_h=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        wind_h(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    wind_h(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,53);
    wind_h(lon_ind(i,1),lat_ind(i,1),12)=2;
end

RGL=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        RGL(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    RGL(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,54);
end

rad_atten=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        rad_atten(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    rad_atten(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,55);
end

wind_atten=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)   
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        wind_atten(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    wind_atten(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,56);
end

trunk_ratio=nan(length(lon),length(lat),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        trunk_ratio(lon_ind(i,1),lat_ind(i,1),:)=0;
    end
    trunk_ratio(lon_ind(i,1),lat_ind(i,1),1:11)=veg_lib(:,57);
end

albedo=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        albedo(lon_ind(i,1),lat_ind(i,1),:,:)=0;
    end
    for j=1:length(veg_class)-1
        albedo(lon_ind(i,1),lat_ind(i,1),:,j)=veg_lib(j,17:28);
    end
    albedo(lon_ind(i,1),lat_ind(i,1),:,12)=0.2;
end

veg_rough=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        veg_rough(lon_ind(i,1),lat_ind(i,1),:,:)=0;
    end
    for j=1:length(veg_class)-1
        veg_rough(lon_ind(i,1),lat_ind(i,1),:,j)=veg_lib(j,29:40);
    end
    veg_rough(lon_ind(i,1),lat_ind(i,1),:,12)=0;
end

displacement=nan(length(lon),length(lat),length(month),length(veg_class));
for i=1:length(lon_ind)
    if ~isnan(mask(lon_ind(i,1),lat_ind(i,1)))
        displacement(lon_ind(i,1),lat_ind(i,1),:,:)=0;
    end
    for j=1:length(veg_class)-1
        displacement(lon_ind(i,1),lat_ind(i,1),:,j)=veg_lib(j,41:52);
    end
    displacement(lon_ind(i,1),lat_ind(i,1),:,12)=0;
end


end
% ##########################################################################


% ##########################################################################      
for jj=1:1
nccreate('Sierra_VIC5_params.nc','lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

      
nccreate('Sierra_VIC5_params.nc','layer',...
    'Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},...
          'Format','classic')
          
nccreate('Sierra_VIC5_params.nc','run_cell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','gridcell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','lats',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','lons',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','infilt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Ds',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Dsmax',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Ws',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','c',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','expt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Ksat',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','phi_s',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','init_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','elev',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','avg_T',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','dp',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')      

nccreate('Sierra_VIC5_params.nc','bubble',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','quartz',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','bulk_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','soil_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','off_gmt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Wcr_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Wpwp_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','snow_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')      

nccreate('Sierra_VIC5_params.nc','annual_prec',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','resid_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','fs_active',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','snow_band',...
    'Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','cellnum',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','AreaFract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','elevation',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')
     
nccreate('Sierra_VIC5_params.nc','Pfactor',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')
     
nccreate('Sierra_VIC5_params.nc','veg_class',...
    'Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','veg_descr',...
    'Datatype','char',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','root_zone',...
    'Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','month',...
    'Datatype','int32',...
    'Dimensions',{'month',length(month)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Nveg',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','Cv',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','root_depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','root_fract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','LAI',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','overstory',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','rarc',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','rmin',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','wind_h',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','RGL',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','rad_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','wind_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','trunk_ratio',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','albedo',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate('Sierra_VIC5_params.nc','veg_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate('Sierra_VIC5_params.nc','displacement',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')   
end      
% ##########################################################################


% ##########################################################################
for jj=1:1
    
ncwrite('Sierra_VIC5_params.nc','lat',lat);
ncwrite('Sierra_VIC5_params.nc','lon',lon);
ncwrite('Sierra_VIC5_params.nc','mask',mask);
ncwrite('Sierra_VIC5_params.nc','layer',layer);
ncwrite('Sierra_VIC5_params.nc','run_cell',run_cell);
ncwrite('Sierra_VIC5_params.nc','gridcell',gridcell);
ncwrite('Sierra_VIC5_params.nc','lats',lats);
ncwrite('Sierra_VIC5_params.nc','lons',lons);
ncwrite('Sierra_VIC5_params.nc','infilt',infilt);
ncwrite('Sierra_VIC5_params.nc','Ds',Ds);
ncwrite('Sierra_VIC5_params.nc','Dsmax',Dsmax);
ncwrite('Sierra_VIC5_params.nc','Ws',Ws);
ncwrite('Sierra_VIC5_params.nc','c',c);
ncwrite('Sierra_VIC5_params.nc','expt',expt);
ncwrite('Sierra_VIC5_params.nc','Ksat',Ksat);
ncwrite('Sierra_VIC5_params.nc','phi_s',phi_s);
ncwrite('Sierra_VIC5_params.nc','init_moist',init_moist);
ncwrite('Sierra_VIC5_params.nc','elev',elev);
ncwrite('Sierra_VIC5_params.nc','depth',depth);
ncwrite('Sierra_VIC5_params.nc','avg_T',avg_T);
ncwrite('Sierra_VIC5_params.nc','dp',dp);
ncwrite('Sierra_VIC5_params.nc','bubble',bubble);
ncwrite('Sierra_VIC5_params.nc','quartz',quartz);
ncwrite('Sierra_VIC5_params.nc','bulk_density',bulk_density);
ncwrite('Sierra_VIC5_params.nc','soil_density',soil_density);
ncwrite('Sierra_VIC5_params.nc','off_gmt',off_gmt);
ncwrite('Sierra_VIC5_params.nc','Wcr_FRACT',Wcr_FRACT);
ncwrite('Sierra_VIC5_params.nc','Wpwp_FRACT',Wpwp_FRACT);
ncwrite('Sierra_VIC5_params.nc','rough',rough);
ncwrite('Sierra_VIC5_params.nc','snow_rough',snow_rough);   
ncwrite('Sierra_VIC5_params.nc','annual_prec',annual_prec);
ncwrite('Sierra_VIC5_params.nc','resid_moist',resid_moist);
ncwrite('Sierra_VIC5_params.nc','fs_active',fs_active);
ncwrite('Sierra_VIC5_params.nc','snow_band',snow_band);
ncwrite('Sierra_VIC5_params.nc','cellnum',cellnum);
ncwrite('Sierra_VIC5_params.nc','AreaFract',AreaFract);
ncwrite('Sierra_VIC5_params.nc','elevation',elevation);
ncwrite('Sierra_VIC5_params.nc','Pfactor',Pfactor);
ncwrite('Sierra_VIC5_params.nc','veg_class',veg_class);
disp('1')
% ncwrite('Sierra_VIC5_params.nc','veg_descr',veg_descr);
disp('2')
ncwrite('Sierra_VIC5_params.nc','root_zone',root_zone);
ncwrite('Sierra_VIC5_params.nc','month',month);
ncwrite('Sierra_VIC5_params.nc','Nveg',Nveg);
ncwrite('Sierra_VIC5_params.nc','Cv',Cv);
ncwrite('Sierra_VIC5_params.nc','root_depth',root_depth);
ncwrite('Sierra_VIC5_params.nc','root_fract',root_fract);
ncwrite('Sierra_VIC5_params.nc','LAI',LAI);
ncwrite('Sierra_VIC5_params.nc','overstory',overstory);
ncwrite('Sierra_VIC5_params.nc','rarc',rarc);
ncwrite('Sierra_VIC5_params.nc','rmin',rmin);
ncwrite('Sierra_VIC5_params.nc','wind_h',wind_h);
ncwrite('Sierra_VIC5_params.nc','RGL',RGL);
ncwrite('Sierra_VIC5_params.nc','rad_atten',rad_atten);
ncwrite('Sierra_VIC5_params.nc','wind_atten',wind_atten);
ncwrite('Sierra_VIC5_params.nc','trunk_ratio',trunk_ratio);
ncwrite('Sierra_VIC5_params.nc','albedo',albedo);
ncwrite('Sierra_VIC5_params.nc','veg_rough',veg_rough);
ncwrite('Sierra_VIC5_params.nc','displacement',displacement);

end
% ##########################################################################      



ncdisp('./shasta_vic5/Shasta_VIC5_params.nc')

a=[];
a=ncread('./shasta_vic5/Shasta_VIC5_params.nc','gridcell');

ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc')
run=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc','run_cell');
gd=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc','gridcell');

ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc')



