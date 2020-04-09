% Subset parameter
%
% Function for subsetting the VIC 5 image mode parameter file
% Warning: can use a lot of RAM.
%
% INPUTS
% Extent = study area extent. Can input either a shapefile or a list of latlon coordinates
% global_domain = name of the input domain file you wish to subset`
% outname = name of the output, subsetted domain file
%
% OUTPUTS
% Subsetted parameter file for the VIC 5 image driver
%
% Sample inputs:
% extent = horzcat([-121.5942; -118.9067], [37.4058; 38.4058]);
% global_params = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_3/VICGlobal_params.nc';
% outname = '/Volumes/HD4/SWOTDA/Data/Tuolumne/params_sub_1.nc';

function outname = subset_parameter(extent, global_params, outname)

% Read coordinates for full domain
lat = ncread(global_params,'lat');
lon = ncread(global_params,'lon');

if ischar(extent)
    tmp1 = strsplit(extent, '.');
    extension = tmp1{2};
    if strcmp(extension, 'shp')
        disp('Input is a shapefile');
        extent = shaperead(extent);
        lon_sub = extent.X(1:end-1)';
        lat_sub = extent.Y(1:end-1)'; 
    else
        disp('Please input extent as a shapefile or a list of coordinates')
    end
end

if isnumeric(extent)
    lon_sub = extent(:,1);
    lat_sub = extent(:,2);
    disp('Input is a list of coordinates');
end

min_lat=min(lat_sub);
max_lat=max(lat_sub);
min_lon=min(lon_sub);
max_lon=max(lon_sub);

% inclusive (larger domain)
resolution = 1/16; 
disp('Assuming resolution is 1/16 degrees')
lat_ind1 = find(abs(lat(:,1) - min_lat) <= resolution/2);
lat_ind2 = find(abs(lat(:,1) - max_lat) <= resolution/2);
lat_ind = lat_ind1:lat_ind2;
lon_ind1 = find(abs(lon(:,1) - min_lon) <= resolution/2);
lon_ind2 = find(abs(lon(:,1) - max_lon) <= resolution/2);
lon_ind = lon_ind1:lon_ind2;

% exclusive (smaller domain)
% lat_ind = find(lat(:,1)>=min_lat & lat(:,1)<=max_lat);
% lon_ind = find(lon(:,1)>=min_lon & lon(:,1)<=max_lon);

if isempty(lat_ind)
    disp('empty lat ind')
    outname = 'none';
    return
end

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

nccreate(outname,'lat','Datatype','double',...
    'Dimensions',{'lat',length(lat_ind)},'Format','netcdf4_classic')
nccreate(outname,'lon','Datatype','double',...
    'Dimensions',{'lon',length(lon_ind)},'Format','netcdf4_classic')
nccreate(outname,'layer','Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},'Format','netcdf4_classic')
nccreate(outname,'snow_band','Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},'Format','netcdf4_classic')
nccreate(outname,'veg_class','Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},'Format','netcdf4_classic')
nccreate(outname,'root_zone','Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},'Format','netcdf4_classic')
nccreate(outname,'month','Datatype','int32',...
    'Dimensions',{'month',length(month)},'Format','netcdf4_classic')

ncwrite(outname,'lat',lat(lat_ind));
ncwrite(outname,'lon',lon(lon_ind));
ncwrite(outname,'layer',layer);
ncwrite(outname,'snow_band',snow_band);
ncwrite(outname,'veg_class',veg_class);
ncwrite(outname,'root_zone',root_zone);
ncwrite(outname,'month',month);

% ##########################################################################
% 2-D matrices
var_list2={'mask','lats','lons','infilt','Ds','Dsmax','Ws','c',...
    'elev','avg_T','dp','off_gmt','rough','snow_rough','annual_prec','July_Tavg',...
    'cellnum'};
for i=1:length(var_list2)
    var_in=[];
    var_in=ncread(global_params,var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate(outname,var_list2{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list2{i},eval(var_list2{i}));
end

% Several variables must be integers, not floats
% These are run_cell, gridcell, and fs_active
% 2-D matrices
var_list2={'run_cell','gridcell','fs_active', 'Nveg'};
for i=1:length(var_list2)
    var_in=[];
    var_in=ncread(global_params,var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate(outname,var_list2{i},...
        'Datatype','int32',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list2{i},eval(var_list2{i}));
end


% ##########################################################################
% 3-D matrices
var_list3={'expt','Ksat','phi_s','init_moist','depth','bubble','quartz','bulk_density',...
    'soil_density','Wcr_FRACT','Wpwp_FRACT','resid_moist','AreaFract','elevation',...
    'Pfactor','Cv','overstory','rarc','rmin','wind_h','RGL','rad_atten','wind_atten','trunk_ratio'};

% integer variables: overstory
var_in = ncread(global_params,'overstory');
overstory = var_in(lon_ind, lat_ind, :);
nccreate(outname,'overstory','Datatype','int32','Dimensions',...
    {'lon',length(lon_ind),'lat',length(lat_ind),'veg_class',length(veg_class)},...
    'Format','netcdf4_classic')   
ncwrite(outname,'overstory',overstory);

for i=1:12
    var_in=[];
    var_in=ncread(global_params,var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate(outname,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'nlayer',length(layer)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list3{i},eval(var_list3{i}));
end
for i=13:15
    var_in=[];
    var_in=ncread(global_params,var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate(outname,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'snow_band',length(snow_band)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list3{i},eval(var_list3{i}));
end
for i=16:24
    
    if i==17 % skip overstory (int, not double)
        continue
    end
    
    var_in=[];
    var_in=ncread(global_params,var_list3{i});
    eval([var_list3{i},'=var_in(lon_ind,lat_ind,:);']);
    
    nccreate(outname,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list3{i},eval(var_list3{i}));
end

% ##########################################################################
% 4-D matrices

var_list4={'root_depth','root_fract','LAI','albedo','veg_rough','displacement','fcanopy'};
for i=1:2
    var_in=[];
    var_in=ncread(global_params,var_list4{i});
    eval([var_list4{i},'=var_in(lon_ind,lat_ind,:,:);']);
    
    nccreate(outname,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list4{i},eval(var_list4{i}));
end
for i=3:6
    var_in=[];
    var_in=ncread(global_params,var_list4{i});
    eval([var_list4{i},'=var_in(lon_ind,lat_ind,:,:);']);
    
    nccreate(outname,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'month',length(month),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list4{i},eval(var_list4{i}));
end
for i=7
    var_in=[];
    var_in=ncread(global_params,'fcanopy');
    eval([var_list4{i},'=var_in(lon_ind,lat_ind,:,:);']);
    
    nccreate(outname,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind),'month',length(month),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list4{i},eval(var_list4{i}));
end

return