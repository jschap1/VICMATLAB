% Modifying the parameter file
%
% LAI and fcanopy values for water are currently set to 0. However, VIC
% cannot handle values of zero; it results in NaN. Therefore, this script
% modifies the image driver's parameter file to set LAI and fcanopy to
% small positive values, instead of 0.

orig_params = '/Volumes/HD_ExFAT/output/VICGlobal_params.nc';
mod_params = '/Volumes/HD_ExFAT/output/VICGlobal_params_mod.nc';

% Read coordinates
lat = ncread(orig_params,'lat');
lon = ncread(orig_params,'lon');

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

nccreate(mod_params,'lat','Datatype','double',...
    'Dimensions',{'lat',length(lat)},'Format','netcdf4_classic')
nccreate(mod_params,'lon','Datatype','double',...
    'Dimensions',{'lon',length(lon)},'Format','netcdf4_classic')
nccreate(mod_params,'layer','Datatype','int32',...
    'Dimensions',{'nlayer',length(layer)},'Format','netcdf4_classic')
nccreate(mod_params,'snow_band','Datatype','int32',...
    'Dimensions',{'snow_band',length(snow_band)},'Format','netcdf4_classic')
nccreate(mod_params,'veg_class','Datatype','int32',...
    'Dimensions',{'veg_class',length(veg_class)},'Format','netcdf4_classic')
nccreate(mod_params,'root_zone','Datatype','int32',...
    'Dimensions',{'root_zone',length(root_zone)},'Format','netcdf4_classic')
nccreate(mod_params,'month','Datatype','int32',...
    'Dimensions',{'month',length(month)},'Format','netcdf4_classic')

ncwrite(mod_params,'lat',lat);
ncwrite(mod_params,'lon',lon);
ncwrite(mod_params,'layer',layer);
ncwrite(mod_params,'snow_band',snow_band);
ncwrite(mod_params,'veg_class',veg_class);
ncwrite(mod_params,'root_zone',root_zone);
ncwrite(mod_params,'month',month);

% ##########################################################################
% 2-D matrices
var_list2={'mask','lats','lons','infilt','Ds','Dsmax','Ws','c',...
    'elev','avg_T','dp','off_gmt','rough','snow_rough','annual_prec','July_Tavg',...
    'cellnum'};
for i=1:length(var_list2)
    var_in=[];
    var_in=ncread(orig_params,var_list2{i});
    eval([var_list2{i},'=var_in;']);
    
    nccreate(mod_params,var_list2{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list2{i},eval(var_list2{i}));
end

% Several variables must be integers, not floats
% These are run_cell, gridcell, and fs_active
% 2-D matrices
var_list2={'run_cell','gridcell','fs_active', 'Nveg'};
for i=1:length(var_list2)
    var_in=[];
    var_in=ncread(orig_params,var_list2{i});
    eval([var_list2{i},'=var_in;']);
    
    nccreate(mod_params,var_list2{i},...
        'Datatype','int32',...
        'Dimensions',{'lon',length(lon),'lat',length(lat)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list2{i},eval(var_list2{i}));
end


% ##########################################################################
% 3-D matrices
var_list3={'expt','Ksat','phi_s','init_moist','depth','bubble','quartz','bulk_density',...
    'soil_density','Wcr_FRACT','Wpwp_FRACT','resid_moist','AreaFract','elevation',...
    'Pfactor','Cv','overstory','rarc','rmin','wind_h','RGL','rad_atten','wind_atten','trunk_ratio'};

% integer variables: overstory
var_in = ncread(orig_params,'overstory');
overstory = var_in;
nccreate(mod_params,'overstory','Datatype','int32','Dimensions',...
    {'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
    'Format','netcdf4_classic')   
ncwrite(mod_params,'overstory',overstory);

for i=1:12
    var_in=[];
    var_in=ncread(orig_params,var_list3{i});
    eval([var_list3{i},'=var_in;']);
    
    nccreate(mod_params,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'nlayer',length(layer)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list3{i},eval(var_list3{i}));
end
for i=13:15
    var_in=[];
    var_in=ncread(orig_params,var_list3{i});
    eval([var_list3{i},'=var_in;']);
    
    nccreate(mod_params,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'snow_band',length(snow_band)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list3{i},eval(var_list3{i}));
end
for i=16:24
    
    if i==17 % skip overstory (int, not double)
        continue
    end
    
    var_in=[];
    var_in=ncread(orig_params,var_list3{i});
    eval([var_list3{i},'=var_in;']);
    
    nccreate(mod_params,var_list3{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list3{i},eval(var_list3{i}));
end

% ##########################################################################
% 4-D matrices

var_list4={'root_depth','root_fract','LAI','albedo','veg_rough','displacement','fcanopy'};
for i=1:2
    var_in=[];
    var_in=ncread(orig_params,var_list4{i});
    eval([var_list4{i},'=var_in;']);
    
    nccreate(mod_params,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'root_zone',length(root_zone),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list4{i},eval(var_list4{i}));
end
for i=3:6
    var_in=[];
    var_in=ncread(orig_params,var_list4{i});
    
    if i==3 % changing the value of LAI
        water_ind = find(var_in(:,:,:,1)==0);
        water_LAI = var_in(:,:,:,1);
        water_LAI(water_ind) = 1e-3;
        var_in(:,:,:,1) = water_LAI;
    end
    
    eval([var_list4{i},'=var_in;']);
    
    nccreate(mod_params,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list4{i},eval(var_list4{i}));
end
for i=7
    var_in=[];
    var_in=ncread(orig_params,'fcanopy');
    
    % change the value of fcanopy
    water_ind = find(var_in(:,:,:,1)==0);
    water_fcan = var_in(:,:,:,1);
    water_fcan(water_ind) = 1e-3;
    var_in(:,:,:,1) = water_fcan;    
    
    eval([var_list4{i},'=var_in;']);
    
    nccreate(mod_params,var_list4{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon),'lat',length(lat),'month',length(month),'veg_class',length(veg_class)},...
        'Format','netcdf4_classic')    
    ncwrite(mod_params,var_list4{i},eval(var_list4{i}));
end
