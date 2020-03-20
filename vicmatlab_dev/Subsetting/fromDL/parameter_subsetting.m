%% domain subsetting

clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IN THIS PART THE READERS NEED TO PUT IN THE DOMAIN EXTENT, EITHER SHAPEFILE OR A LIST OF
% LAT/LON. THEN THE PROGRAM CAN FIGURE OUT LAT_IND & LON_IND

% min_lat=min(lat_rec);
% max_lat=max(lat_rec);
% min_lon=min(lon_rec);
% max_lon=max(lon_rec);

ori_dir='/Volumes/Disk4/VICGlobal/output/';
out_dir='/Volumes/Disk4/VICGlobal/output/';
lat=ncread([ori_dir,'VICGlobal_params.nc'],'lat');
lon=ncread([ori_dir,'VICGlobal_params.nc'],'lon');

lat_ind=find(lat(:,1)>=37.6 & lat(:,1)<=38.3);
lon_ind=find(lon(:,1)>=-120.6 & lon(:,1)<=-119.2);
% ##########################################################################

% 1-D array
var_list1={'lat','lon'};
nccreate([out_dir,'subset_VICGlobal_domain.nc'],'lat','Datatype','double',...
    'Dimensions',{'lat',length(lat_ind)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_domain.nc'],'lon','Datatype','double',...
    'Dimensions',{'lon',length(lon_ind)},'Format','netcdf4_classic')
ncwrite([out_dir,'subset_VICGlobal_domain.nc'],'lat',lat(lat_ind));
ncwrite([out_dir,'subset_VICGlobal_domain.nc'],'lon',lon(lon_ind));

% 2-D matrices
var_list2={'mask','frac','area'};

i=1;
var_in=[];
var_in=ncread([ori_dir,'VICGlobal_domain.nc'],var_list2{i});
eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);

nccreate([out_dir,'subset_VICGlobal_domain.nc'],var_list2{i},...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
    'Format','netcdf4_classic')    
ncwrite([out_dir,'subset_VICGlobal_domain.nc'],var_list2{i},eval(var_list2{i}));

for i=2:length(var_list2)
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_domain.nc'],var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate([out_dir,'subset_VICGlobal_domain.nc'],var_list2{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_domain.nc'],var_list2{i},eval(var_list2{i}));
end

%% Parameter subsetting

clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IN THIS PART THE READERS NEED TO PUT IN THE DOMAIN EXTENT, EITHER SHAPEFILE OR A LIST OF
% LAT/LON. THEN THE PROGRAM CAN FIGURE OUT LAT_IND & LON_IND

% min_lat=min(lat_rec);
% max_lat=max(lat_rec);
% min_lon=min(lon_rec);
% max_lon=max(lon_rec);

ori_dir='/Volumes/Disk4/VICGlobal/output/';
out_dir='/Volumes/Disk4/VICGlobal/output/';
lat=ncread([ori_dir,'VICGlobal_params.nc'],'lat');
lon=ncread([ori_dir,'VICGlobal_params.nc'],'lon');

lat_ind=find(lat(:,1)>=37.6 & lat(:,1)<=38.3);
lon_ind=find(lon(:,1)>=-120.6 & lon(:,1)<=-119.2);

% ##########################################################################
% 1-D array 
var_list1={'lat','lon','layer','snow_band','veg_class','veg_descr','root_zone','month'};

layer=[1;2;3];
snow_band=[1;2;3;4;5];
veg_class=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17];
veg_descr=char({'Open Water','Evergreen Needleleaf','Evergreen Broadleaf','Deciduous Needleleaf',...
    'Deciduous Broadleaf','Mixed Forest','Closed Shrub','Open Shrub','Savannas','Woody Savannas',...
    'Grasslands','Permanent Wetlands','Crop land','Urban','Cropland NAtural Veg Mosaic',...
    'SnowIce','Barren'});
root_zone=[1;2;3];
month=[1;2;3;4;5;6;7;8;9;10;11;12];

nccreate([out_dir,'subset_VICGlobal_params.nc'],'lat','Datatype','double',...
    'Dimensions',{'lat',length(lat_ind)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'lon','Datatype','double',...
    'Dimensions',{'lon',length(lon_ind)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'layer','Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'snow_band','Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'veg_class','Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'root_zone','Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},'Format','netcdf4_classic')
nccreate([out_dir,'subset_VICGlobal_params.nc'],'month','Datatype','int32',...
    'Dimensions',{'month',length(month)},'Format','netcdf4_classic')

ncwrite([out_dir,'subset_VICGlobal_params.nc'],'lat',lat(lat_ind));
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'lon',lon(lon_ind));
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'layer',layer);
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'snow_band',snow_band);
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'veg_class',veg_class);
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'root_zone',root_zone);
ncwrite([out_dir,'subset_VICGlobal_params.nc'],'month',month);

% ##########################################################################
% 2-D matrices
var_list2={'mask','run_cell','gridcell','fs_active','Nveg','lats','lons','infilt','Ds','Dsmax','Ws','c',...
    'elev','avg_T','dp','off_gmt','rough','snow_rough','annual_prec','July_Tavg',...
    'cellnum'};

for i=1:5
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list2{i},...
        'Datatype','int32',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list2{i},eval(var_list2{i}));
end

for i=6:length(var_list2)
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list2{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list2{i},eval(var_list2{i}));
end

% ##########################################################################
% 3-D matrices
var_list3={'expt','Ksat','phi_s','init_moist','depth','bubble','quartz','bulk_density',...
    'soil_density','Wcr_FRACT','Wpwp_FRACT','resid_moist','AreaFract','elevation',...
    'Pfactor','overstory','Cv','rarc','rmin','wind_h','RGL','rad_atten','wind_atten','trunk_ratio'};

for i=1:12
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'nlayer',length(layer)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},eval(var_list3{i}));
end
for i=13:15
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'snow_band',length(snow_band)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},eval(var_list3{i}));
end

i=16;
var_in=[];
var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list3{i});
eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'veg_class',length(veg_class)},...
    'Format','netcdf4_classic')    
ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},eval(var_list3{i}));

for i=17:24
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list3{i},eval(var_list3{i}));
end

% ##########################################################################
% 4-D matrices

var_list4={'root_depth','root_fract','LAI','albedo','veg_rough','displacement','fcanopy'};
for i=1:2
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list4{i});
    eval([var_list4{i},'=var_in(lon_ind,lat_ind,:,:);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list4{i},eval(var_list4{i}));
end
for i=3:7
    var_in=[];
    var_in=ncread([ori_dir,'VICGlobal_params.nc'],var_list4{i});
    eval([var_list4{i},'=var_in(lon_ind,lat_ind,:,:);']);
    
    nccreate([out_dir,'subset_VICGlobal_params.nc'],var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'month',length(month),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite([out_dir,'subset_VICGlobal_params.nc'],var_list4{i},eval(var_list4{i}));
end




























