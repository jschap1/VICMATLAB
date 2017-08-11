% Generic workflow for creating VIC input files
%
% Clips daily forcing and soil parameter data for CONUS to just the cells within 
% a specified basin shapefile. Saves the subsetted forcing and soil
% parameter data in an appropriate format to use for VIC input.
%
% Is set up to use Livneh forcing data, which is on a 444 by 922 lat/lon
% grid (1/16 degree resolution)

%% Specify inputs

% Directory of daily CONUS met. forcing file
forcingpath = '/Users/jschapMac/Documents/HydrologyData/Livneh/MetNC';

% Directory where clipped forcing files should be saved
forcingsavedir = '/Users/jschapMac/Desktop/Tuolumne2/Forcings';

% Number of forcings in the daily CONUS met. forcing file
numforcings = 4;

% Beginning and ending years of simulation (must be included in the daily CONUS
% met. forcing file)
beginyear = 2006;
endyear = 2011;

% Number of decimal points of precision to use for forcing file names
grid_decimal = 5;

% Directory of CONUS soil parameter file
soilpath = '/Users/jschapMac/Documents/HydrologyData/VICParametersCONUS';
soilname = 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt'; 

% Directory where clipped soil parameter file should be saved
soilsavedir = '/Users/jschapMac/Desktop/Tuolumne2';

%% Load the mask

addpath(forcingpath)
metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');

%%%
% Run this code to convert lon coords if they use E/W, 
% instead of absolute value system
metlon = metlon - 360;
%%%

% Load coarse resolution basin mask (1/16 deg., same as the forcing data)
[mask, R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne/RoutingInputs/basinmask_coarse.asc');
[ind1, ind2] = find(mask);

% Get lat/lon of basin mask (only the pixels whose values are 1)
ncells = sum(mask(:)~=0);
masklat = NaN(ncells,1);
masklon = NaN(ncells,1);

for i=1:ncells
    [masklat(i),masklon(i)] = pix2latlon(R, ind1(i), ind2(i));
end

%% Extract the forcings whose lat/lon match those of the basin mask

disp('Extracting forcing variables')

nyears = endyear - beginyear + 1;
t_ind = 1;
cum_days = 0;

for t = beginyear:endyear
    
    if t==beginyear, tic, end
    prec = ncread(['prec.' num2str(t) '.nc'], 'prec');
    tmax = ncread(['tmax.' num2str(t) '.nc'], 'tmax');
    tmin = ncread(['tmin.' num2str(t) '.nc'], 'tmin');
    wind = ncread(['wind.' num2str(t) '.nc'], 'wind');
       
    info = ncinfo(['prec.' num2str(t) '.nc']);
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells);
    
    for k=1:ncells
        % Get the index of the met. forcing data that matches the basin mask
        [Lia,lat_ind] = ismember(masklat(k),metlat);
        [Lia,lon_ind] = ismember(masklon(k),metlon);
        data(:,1,k) = prec(lon_ind,lat_ind, :);        
        data(:,2,k) = tmin(lon_ind,lat_ind, :);
        data(:,3,k) = tmax(lon_ind,lat_ind, :);
        data(:,4,k) = wind(lon_ind,lat_ind, :);
    end
    
    if t_ind~=1
        data_cum = vertcat(data_cum, data);
    else
        data_cum = data;
    end

    t_ind = t_ind + 1;
    cum_days = size(prec,3) + cum_days;
    
    if t==beginyear 
        disp(['About ' num2str(toc*nyears/60) ...
            ' minutes remaining.'])
    end
    
end

% % Save met. forcings as .mat file
% % The dimensions are [numdays, numforcings, ncells]
% save('METFORC.mat', 'data_cum');

fstring = ['%.' num2str(grid_decimal) 'f'];
for k=1:ncells     
    filename = ['data_' num2str(masklat(k),fstring) '_' num2str(masklon(k),fstring)];
    dlmwrite(fullfile(forcingsavedir, filename), data_cum(:,:,i))  
end
display(['Forcing data saved to ' forcingsavedir])

%% Extract soils data

% Load the soils data

soils = load(fullfile(soilpath, soilname));

slat = soils(:,3);
slon = soils(:,4);

soils(:,1) = 0;

for k=1:ncells
    % Get the index of the met. forcing data that matches the basin mask
    sind = find(masklat(k) == slat & masklon(k) == slon);
    soils(sind,1) = 1;
end

fstring = ['%.' num2str(grid_decimal) 'f'];
fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];
fID = fopen(fullfile(soilsavedir, 'soilsl.TUOLUMNE'),'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved to ' soilsavedir])

% Note: the delimiter and the format spec must be specified precisely as
% the Stehekin example from the VIC website in order to avoid the error
% about CELL LATITUDE not being found/for VIC to successfully read the soil
% parameter file.

%% Make plots to check that everything is all right

% Plot basin mask
figure 
mapshow(mask,R,'DisplayType','surface')

% Plot met. forcing file lat lons
% figure
% plot(data_cum(:,:,1),'.')

% Plot soils file lat lons
figure
plot(soilsclip(:,4), soilsclip(:,3), '.')

% Plot basin mask and soil lat lons in one figure:
% figure,hold on, mapshow(mask,R),  mapshow(soilsclip(:,4), soilsclip(:,3))