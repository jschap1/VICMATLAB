% Converts vegetation parameter files to geotiff files
%
% INPUTS
% Outputs from load_veg_parameters
% veglib = name of vegetation library file. Must be in numeric only format. Can leave blank.
%
% OUTPUTS
% geotiffs containing:
% -cover fraction (1*ntypes) for each vegetation type
% -rooting depth 1 and 2 (2*ntypes)
% -rooting fraction 1 and 2 (2*ntypes)
% -LAI for each month (12*ntypes)
% -nveg (1*ntypes)
% -landmask (1*ntypes)

function [] = convert_veg_parameters_v2(nvegtable, vegparamtable, latlontable, veglib, outdir)

lat = latlontable(:,2);
lon = latlontable(:,3);

nlat = length(unique(lat));
nlon = length(unique(lon));

rasterSize = [nlat, nlon];
latlim = [min(lat), max(lat)];
lonlim = [min(lon), max(lon)];

R = georefcells(latlim,lonlim,rasterSize);

%% Create nveg map and landmask

ncells = length(lat);
mask_vect = nvegtable(:,2)>0;
nveg_map = xyz2grid(lon, lat, nvegtable(:,2));
vegmask = xyz2grid(lon, lat, mask_vect);
geotiffwrite(fullfile(outdir, 'nveg_map.tif'), flipud(nveg_map), R)
geotiffwrite(fullfile(outdir, 'veg_mask.tif'), flipud(vegmask), R)

figure, imagesc(nveg_map), title('number of vegetation classes')
figure, imagesc(vegmask), title('vegetation mask')

%% Create cover fraction maps for each vegetation type

vegtype = vegparamtable(:,2); % vegetation type
vegfract = vegparamtable(:,3); % cover fraction
rootfract = vegparamtable(:,4:5); % root fraction
rootdepth = vegparamtable(:,6:7); % root depth
LAI = vegparamtable(:,8:19); % monthly LAI

[~, ind] = ismember(vegparamtable(:,1), nvegtable(:,1));

vegtype_vect = zeros(ncells, 1);
vegfract_vect = zeros(ncells, 1);
rootfract_arr = zeros(ncells, 2);
rootdepth_arr = zeros(ncells, 2);
LAI_arr = zeros(ncells, 12);

vegtype_vect(ind) = vegtype;
vegfract_vect(ind) = vegfract;
rootfract_arr(ind, :) = rootfract;
rootdepth_arr(ind, :) = rootdepth;
LAI_arr(ind, :) = LAI;

vegtype_map = xyz2grid(lon, lat, vegtype_vect);
vegfract_map = xyz2grid(lon, lat, vegfract_vect);

rootfract1_map = xyz2grid(lon, lat, rootfract_arr(:,1));
rootfract2_map = xyz2grid(lon, lat, rootfract_arr(:,2));

rootdepth1_map = xyz2grid(lon, lat, rootdepth_arr(:,1));
rootdepth2_map = xyz2grid(lon, lat, rootdepth_arr(:,2));

LAI1_map = xyz2grid(lon, lat, LAI_arr(:,1));

figure
subplot(3,2,1), imagesc(vegtype_map), title('vegetation types')
subplot(3,2,2), imagesc(vegfract_map), title('cover fraction')
subplot(3,2,3), imagesc(rootfract1_map), title('root fraction 1')
subplot(3,2,4), imagesc(rootdepth1_map), title('root depth 1')
subplot(3,2,5), imagesc(LAI1_map), title('January LAI')

% END

%% Use vegetation library to make maps of other variables

if ~isempty(veglib)
    
    vegetation_library = load(veglib);
    
    % albedo
    mean_albedo_per_vegtype = mean(vegetation_library(:,17:28), 2);
    
    vegmask = xyz2grid(lon, lat, mask_vect);
    geotiffwrite(fullfile(outdir, 'nveg_map.tif'), flipud(vegmask), R)
    
    % calculate albedo for the vegetation type
    albedo_vect = zeros(size(vegtype_vect));
    nvegtypes = 11;
    for i=1:nvegtypes
        albedo_vect(vegtype_vect == i) = mean_albedo_per_vegtype(i);
    end
    
    albedo_map = xyz2grid(lon, lat, albedo_vect);
    outname = fullfile(outdir, 'albedo_map.tif');    
    geotiffwrite(outname, flipud(albedo_map), R)

% Need to account for bare soils somehow.
% Also, the albedo values seem kind of low.
    
end

%% Write out to geotiffs

% fname_map = [vegnames{k} '_' months{m} '_LAI_map.tif'];
% geotiffwrite(fname_map, flipud(lai_map), R)

return

