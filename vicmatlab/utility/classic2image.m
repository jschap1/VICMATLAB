% Classic to Image
%
% Converts VIC Classic Driver inputs to VIC Image Driver inputs
% Written by Dongyue Li
% Adapted 4/9/2020 JRS
%
% INPUTS (ASCII text files)
% Forcing files 
% Soil parameter file
% Vegetation parameter file
% Vegetation library file
% Elevation bands file
%
% OUTPUTS (NetCDF)
% Domain file
% Parameter file
% Forcing files
%
% Sample inputs
% inputs.forcdir = 'forcing_474/data*';
% inputs.veglib = '/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/vegetation_library.mat';
% inputs.soilparfile = '';
% inputs.snowband = '';
% inputs.vegparam = '';
% inputs.domainfile_name = 'Shasta_VIC5_domain.nc';
% inputs.params_name = 'Sierra_VIC5_params.nc';

% Notes: 
% This function is hard-coded for a 1/16 degree resolution domain
% This function is hard-coded to work only for the VIC inputs from Livneh
% et al.
% Tested with Upper Tuolumne basin/Livneh parameters 4/9/2020

function outputs = classic2image(inputs)

%% Create domain file

% ##########################################################################
file=dir(inputs.forcdir); % this is the forcing data directory
% file=dir('forcing_474/data*'); % this is the forcing data directory

lat_rec = zeros(length(file), 1);
lon_rec = zeros(length(file), 1);

for i=1:length(file)
    c = strsplit(file(i).name, '_');
    lat_rec(i,1) = str2double(c{end-1});
    lon_rec(i,1) = str2double(c{end});  
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


nccreate(inputs.domainfile_name,'mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')


nccreate(inputs.domainfile_name,'lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','classic')      

nccreate(inputs.domainfile_name,'lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','classic') 
     
nccreate(inputs.domainfile_name,'frac',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate(inputs.domainfile_name,'area',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')


ncwrite(inputs.domainfile_name,'mask',mask);
ncwrite(inputs.domainfile_name,'lon',lon);
ncwrite(inputs.domainfile_name,'lat',lat);
ncwrite(inputs.domainfile_name,'frac',frac);
ncwrite(inputs.domainfile_name,'area',area);

ncdisp(inputs.domainfile_name)

disp('Created domain file. Now creating parameter file.')


% a=[];
% a=ncread(inputs.domainfile_name,'lon');

% % working Sierra Nevada version
% ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc')
% 
% a=[];
% a=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_domain.nc','mask');


%% Create parameter file for the domain (same extent with domain file)

% clear
% clc

% FillValue_int=-2147483647;
% FillValue_f=9.969209968386869e+36;

% ##########################################################################
veg_lib = load(inputs.veglib);
soil=dlmread(inputs.soilparfile);
soil=soil(1:208326,:); % from 208327 to the end of the file is in Canada
snow=dlmread(inputs.snowband);
snow=snow(1:208326,:);
veg_param=dlmread(inputs.vegparam);
% ##########################################################################

% ##########################################################################
% file=dir(inputs.forcdir); % this is the forcing data directory
% lat_rec = zeros(length(file), 1);
% lon_rec = zeros(length(file), 1);
% for i=1:length(file)
%     lat_rec(i,1)=str2double(file(i).name(6:13));    
%     lon_rec(i,1)=str2double(file(i).name(15:24));   
% end

% min_lat=min(lat_rec); % domain boundary
% max_lat=max(lat_rec);
% min_lon=min(lon_rec);
% max_lon=max(lon_rec);
% lat=([min_lat:0.0625:max_lat])'; % domain lat/lon
% lon=([min_lon:0.0625:max_lon])';
% clear min* max*

for i=1:length(file)
    lon_ind(i,1)=find(lon==lon_rec(i,1)); % position of running pix in DOMAIN
    lat_ind(i,1)=find(lat==lat_rec(i,1));
end

soil_ind = zeros(length(file), 1);
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


end % jj (line 181)
% ##########################################################################


% ##########################################################################      
for jj=1:1
nccreate(inputs.params_name,'lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','classic')

nccreate(inputs.params_name,'mask',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

      
nccreate(inputs.params_name,'layer',...
    'Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},...
          'Format','classic')
          
nccreate(inputs.params_name,'run_cell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate(inputs.params_name,'gridcell',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'lats',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'lons',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'infilt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'Ds',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'Dsmax',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'Ws',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'c',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'expt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'Ksat',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')
      
nccreate(inputs.params_name,'phi_s',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'init_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'elev',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'avg_T',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'dp',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')      

nccreate(inputs.params_name,'bubble',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'quartz',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'bulk_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'soil_density',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'off_gmt',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'Wcr_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'Wpwp_FRACT',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate(inputs.params_name,'snow_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')      

nccreate(inputs.params_name,'annual_prec',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'resid_moist',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
          'Format','classic')

nccreate(inputs.params_name,'fs_active',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'snow_band',...
    'Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},...
          'Format','classic')

nccreate(inputs.params_name,'cellnum',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')
      
nccreate(inputs.params_name,'AreaFract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')

nccreate(inputs.params_name,'elevation',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')
     
nccreate(inputs.params_name,'Pfactor',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
          'Format','classic')
     
nccreate(inputs.params_name,'veg_class',...
    'Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate(inputs.params_name,'veg_descr',...
    'Datatype','char',...
    'Dimensions',{'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'root_zone',...
    'Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},...
          'Format','classic')
      
nccreate(inputs.params_name,'month',...
    'Datatype','int32',...
    'Dimensions',{'month',length(month)},...
          'Format','classic')

nccreate(inputs.params_name,'Nveg',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','classic')

nccreate(inputs.params_name,'Cv',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'root_depth',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'root_fract',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'LAI',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'overstory',...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'rarc',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'rmin',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'wind_h',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'RGL',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'rad_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate(inputs.params_name,'wind_atten',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'trunk_ratio',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'albedo',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')

nccreate(inputs.params_name,'veg_rough',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')
      
nccreate(inputs.params_name,'displacement',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
          'Format','classic')   
end %  jj (line 565)      
% ##########################################################################


% ##########################################################################
for jj=1:1
    
ncwrite(inputs.params_name,'lat',lat);
ncwrite(inputs.params_name,'lon',lon);
ncwrite(inputs.params_name,'mask',mask);
ncwrite(inputs.params_name,'layer',layer);
ncwrite(inputs.params_name,'run_cell',run_cell);
ncwrite(inputs.params_name,'gridcell',gridcell);
ncwrite(inputs.params_name,'lats',lats);
ncwrite(inputs.params_name,'lons',lons);
ncwrite(inputs.params_name,'infilt',infilt);
ncwrite(inputs.params_name,'Ds',Ds);
ncwrite(inputs.params_name,'Dsmax',Dsmax);
ncwrite(inputs.params_name,'Ws',Ws);
ncwrite(inputs.params_name,'c',c);
ncwrite(inputs.params_name,'expt',expt);
ncwrite(inputs.params_name,'Ksat',Ksat);
ncwrite(inputs.params_name,'phi_s',phi_s);
ncwrite(inputs.params_name,'init_moist',init_moist);
ncwrite(inputs.params_name,'elev',elev);
ncwrite(inputs.params_name,'depth',depth);
ncwrite(inputs.params_name,'avg_T',avg_T);
ncwrite(inputs.params_name,'dp',dp);
ncwrite(inputs.params_name,'bubble',bubble);
ncwrite(inputs.params_name,'quartz',quartz);
ncwrite(inputs.params_name,'bulk_density',bulk_density);
ncwrite(inputs.params_name,'soil_density',soil_density);
ncwrite(inputs.params_name,'off_gmt',off_gmt);
ncwrite(inputs.params_name,'Wcr_FRACT',Wcr_FRACT);
ncwrite(inputs.params_name,'Wpwp_FRACT',Wpwp_FRACT);
ncwrite(inputs.params_name,'rough',rough);
ncwrite(inputs.params_name,'snow_rough',snow_rough);   
ncwrite(inputs.params_name,'annual_prec',annual_prec);
ncwrite(inputs.params_name,'resid_moist',resid_moist);
ncwrite(inputs.params_name,'fs_active',fs_active);
ncwrite(inputs.params_name,'snow_band',snow_band);
ncwrite(inputs.params_name,'cellnum',cellnum);
ncwrite(inputs.params_name,'AreaFract',AreaFract);
ncwrite(inputs.params_name,'elevation',elevation);
ncwrite(inputs.params_name,'Pfactor',Pfactor);
ncwrite(inputs.params_name,'veg_class',veg_class);
disp('1')
% ncwrite(inputs.params_name,'veg_descr',veg_descr);
disp('2')
ncwrite(inputs.params_name,'root_zone',root_zone);
ncwrite(inputs.params_name,'month',month);
ncwrite(inputs.params_name,'Nveg',Nveg);
ncwrite(inputs.params_name,'Cv',Cv);
ncwrite(inputs.params_name,'root_depth',root_depth);
ncwrite(inputs.params_name,'root_fract',root_fract);
ncwrite(inputs.params_name,'LAI',LAI);
ncwrite(inputs.params_name,'overstory',overstory);
ncwrite(inputs.params_name,'rarc',rarc);
ncwrite(inputs.params_name,'rmin',rmin);
ncwrite(inputs.params_name,'wind_h',wind_h);
ncwrite(inputs.params_name,'RGL',RGL);
ncwrite(inputs.params_name,'rad_atten',rad_atten);
ncwrite(inputs.params_name,'wind_atten',wind_atten);
ncwrite(inputs.params_name,'trunk_ratio',trunk_ratio);
ncwrite(inputs.params_name,'albedo',albedo);
ncwrite(inputs.params_name,'veg_rough',veg_rough);
ncwrite(inputs.params_name,'displacement',displacement);

end
% ##########################################################################      

ncdisp(inputs.params_name)

% a=[];
% a=ncread('./shasta_vic5/Shasta_VIC5_params.nc','gridcell');
% 
% ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc')
% run=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc','run_cell');
% gd=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/parameters/Sierra_VIC5_params.nc','gridcell');
% 
% ncdisp('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc')

outputs = 1;

return