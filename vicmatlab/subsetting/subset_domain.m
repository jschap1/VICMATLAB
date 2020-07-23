% Subset domain
%
% Function for subsetting the VIC 5 image mode domain file
%
% INPUTS
% Extent = study area extent. Can input either a shapefile or a list of latlon coordinates
% global_domain = name of the input domain file you wish to subset`
% outname = name of the output, subsetted domain file
%
% OUTPUTS
% Subsetted domain file for the VIC 5 image driver
%
% Sample inputs:
% extent = horzcat([-121.5942; -118.9067], [37.40583; 38.40583]);
% global_domain = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_3/VICGlobal_domain.nc';
% outname = '/Volumes/HD4/SWOTDA/Data/Tuolumne/domain_sub.nc';

function outname = subset_domain(basinmaskname, global_domain, outname)

extent = basin_mask2coordinate_list(basinmaskname);
basinmask = geotiffread2(basinmaskname);

% Read coordinates for full domain
lat = ncread(global_domain,'lat');
lon = ncread(global_domain,'lon');

% if ischar(extent)
%     tmp1 = strsplit(extent, '.');
%     extension = tmp1{2};
%     if strcmp(extension, 'shp')
%         disp('Input is a shapefile');
%         extent = shaperead(extent);
%         lon_sub = extent.X(1:end-1)';
%         lat_sub = extent.Y(1:end-1)';
%         
% %         resolution = 1/16;
% %         disp('Assuming resolution is 1/16 degrees')
% %         lon_sub = min(lon_sub):resolution:max(lon_sub);
% %         lat_sub = min(lat_sub):resolution:max(lat_sub);
%         
%     else
%         disp('Please input extent as a shapefile or a list of coordinates')
%     end
% end

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
% % exclusive (smaller domain)
inclusive = 0; 
if inclusive
    disp('Using inclusive option. Both options seem to work.')
    resolution = 1/16; 
    disp('Assuming resolution is 1/16 degrees')
    lat_ind1 = find(abs(lat(:,1) - min_lat) <= resolution/2);
    lat_ind2 = find(abs(lat(:,1) - max_lat) <= resolution/2);
    lat_ind = lat_ind1:lat_ind2;
    lon_ind1 = find(abs(lon(:,1) - min_lon) <= resolution/2);
    lon_ind2 = find(abs(lon(:,1) - max_lon) <= resolution/2);
    lon_ind = lon_ind1:lon_ind2;
else
    % Need to have lat/lon values on the same grid
    lat_ind = find(lat(:,1)>=min_lat & lat(:,1)<=max_lat);
    lon_ind = find(lon(:,1)>=min_lon & lon(:,1)<=max_lon);
end

if isempty(lat_ind)
    disp('empty lat ind')
    outname = 'none';
    return
end

% ##########################################################################

% 1-D array
% two possible conventions: [lons, lats] or [lon, lat]
lats = 0; % flag for name of variable
if lats
    var_list1={'lats','lons'};
else
    var_list1={'lat','lon'};
end
nccreate(outname,'lat','Datatype','double',...
    'Dimensions',{'lat',length(lat_ind)},'Format','netcdf4_classic')
nccreate(outname,'lon','Datatype','double',...
    'Dimensions',{'lon',length(lon_ind)},'Format','netcdf4_classic')
ncwrite(outname,'lat',lat(lat_ind));
ncwrite(outname,'lon',lon(lon_ind));

% 2-D matrices
var_list_3 = {'mask'};
var_in = [];
var_in=ncread(global_domain,var_list_3{1});

% var_in
% basinmask
% temp = var_in(lon_ind,lat_ind);
% var_list_3{1} = basinmask';

% eval([var_list_3{1},'=var_in(lon_ind,lat_ind);']);

nccreate(outname,var_list_3{1},...
    'Datatype','int32',...
    'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
    'Format','netcdf4_classic')   
ncwrite(outname,var_list_3{1}, ones(size(basinmask')));
% ncwrite(outname,var_list_3{1}, basinmask');
    
var_list2={'frac','area'};
for i=1:length(var_list2)
    var_in=[];
    var_in=ncread(global_domain,var_list2{i});
    eval([var_list2{i},'=var_in(lon_ind,lat_ind);']);
    
    nccreate(outname,var_list2{i},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon_ind),'lat',length(lat_ind)},...
        'Format','netcdf4_classic')    
    ncwrite(outname,var_list2{i},eval(var_list2{i}));
end

return
