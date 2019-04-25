% Converts vegetation parameter files to geotiff files
% 
% Not the most efficient implementation (lots of for loops), but it gets
% the job done as long as the data are not too large. Requires mapping toolbox
%
% INPUTS
% soil parameter file
% processed vegetation parameter file from load_veg_parameters.m
%
% OUTPUTS
% geotiffs containing:
% -cover fraction (1*ntypes) for each vegetation type
% -rooting depth 1 and 2 (2*ntypes)
% -rooting fraction 1 and 2 (2*ntypes)
% -LAI for each month (12*ntypes)
% -nveg (1*ntypes)
% -landmask (1*ntypes)

function [] = convert_veg_parameters(VEGPAR, soils)

soils = soils(:,1:5);
lat = soils(:,3);
lon = soils(:,4);

nlat = length(unique(lat));
nlon = length(unique(lon));

rasterSize = [nlat, nlon];
latlim = [min(lat), max(lat)];
lonlim = [min(lon), max(lon)];

R = georefcells(latlim,lonlim,rasterSize);

%%

% create nveg map and landmask

ncells = length(VEGPAR);
nveg_vect = zeros(ncells,1);

for i=1:ncells
    nveg_vect(i) = VEGPAR(i).nveg; 
end

mask_vect = nveg_vect>0;
nveg_map = xyz2grid(lon, lat, nveg_vect);
vegmask = xyz2grid(lon, lat, mask_vect);
geotiffwrite('nveg_map.tif', flipud(nveg_map), R)
geotiffwrite('veg_mask.tif', flipud(vegmask), R)

figure, imagesc(vegmask), title('vegetation mask')
figure, imagesc(nveg_map), title('number of vegetation classes')

%% Create cover fraction maps for each vegetation type

ntypes = length(fieldnames(VEGPAR))-2; % number of vegetation types in the vegetation library
vegnames = fieldnames(VEGPAR);
vegnames = vegnames(3:end);

cv_mat = zeros(ncells, ntypes);
for i=1:ncells
    for k=1:ntypes

        cv = VEGPAR(i).(vegnames{k}).cv;
        if isempty(cv)
            cv_mat(i,k) = NaN;
        else
            cv_mat(i,k) = cv;
        end

    end
end

for k=1:ntypes
    cv_map = xyz2grid(lon, lat, cv_mat(:,k));
    fname_map = [vegnames{k} '_coverfract_map.tif'];
    geotiffwrite(fname_map, flipud(cv_map), R)
end

clear cv_map cv_mat

%% Create root depth and root fraction maps for each vegetation type

rd1_mat = zeros(ncells, ntypes);
rd2_mat = zeros(ncells, ntypes);
rf1_mat = zeros(ncells, ntypes);
rf2_mat = zeros(ncells, ntypes);

for i=1:ncells
    for k=1:ntypes

        rd = VEGPAR(i).(vegnames{k}).rootdepth;
        if isempty(rd)
            rd1_mat(i,k) = NaN;
            rd2_mat(i,k) = NaN;
        else
            rd1_mat(i,k) = rd(1);
            rd2_mat(i,k) = rd(2);
        end
        
        rf = VEGPAR(i).(vegnames{k}).rootfract;
        if isempty(rf)
            rf1_mat(i,k) = NaN;
            rf2_mat(i,k) = NaN;
        else
            rf1_mat(i,k) = rf(1);
            rf2_mat(i,k) = rf(2);
        end

    end
end

for k=1:ntypes
    rd1_map = xyz2grid(lon, lat, rd1_mat(:,k));
    fname_map = [vegnames{k} '_rootdepth1_map.tif'];
    geotiffwrite(fname_map, flipud(rd1_map), R)
end

for k=1:ntypes
    rd2_map = xyz2grid(lon, lat, rd2_mat(:,k));
    fname_map = [vegnames{k} '_rootdepth2_map.tif'];
    geotiffwrite(fname_map, flipud(rd2_map), R)
end

for k=1:ntypes
    rf1_map = xyz2grid(lon, lat, rf1_mat(:,k));
    fname_map = [vegnames{k} '_rootfract1_map.tif'];
    geotiffwrite(fname_map, flipud(rf1_map), R)
end

for k=1:ntypes
    rf2_map = xyz2grid(lon, lat, rf2_mat(:,k));
    fname_map = [vegnames{k} '_rootfract2_map.tif'];
    geotiffwrite(fname_map, flipud(rf2_map), R)
end

%% Create monthly LAI maps for each vegetation type

LAI = struct();
months = {'jan','feb','mar','apr','may','jun','jul','aug', ...
    'sep','oct','nov','dec'};

for m=1:12
    LAI.(months{m}) = zeros(ncells, ntypes);
end

for m=1:12
    for i=1:ncells
        for k=1:ntypes

            lai = VEGPAR(i).(vegnames{k}).LAI;
            if isempty(lai)
                LAI.(months{m})(i,k) = NaN;
            else
                LAI.(months{m})(i,k) = lai(m);
            end

        end
    end
end

for m=1:12
    for k=1:ntypes
        lai_map = xyz2grid(lon, lat, LAI.(months{m})(:,k));
        fname_map = [vegnames{k} '_' months{m} '_LAI_map.tif'];
        geotiffwrite(fname_map, flipud(lai_map), R)
    end
end

% % for visualizing results
% xvals = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(1));
% yvals = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(2));
% figure, imagesc(xvals, yvals, vegmask)