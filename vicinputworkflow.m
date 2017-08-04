% Generic workflow for creating VIC input files
%
% Clips daily forcing and soil parameter data for CONUS to just the cells within 
% a specified basin shapefile. Saves the subsetted forcing and soil
% parameter data in an appropriate format to use for VIC input.
%
% Is set up to use Livneh forcing data, which is on a 444 by 922 lat/lon grid

%% Specify inputs

% Directory of CONUS soil parameter file
soilpath = '/Users/jschapMac/Documents/HydrologyData/VICParametersCONUS';
soilname = 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt'; 

% Directory where clipped soil parameter file should be saved
soilsavedir = '/Users/jschapMac/Desktop/Tuolumne';

% Directory of daily CONUS met. forcing file
forcingpath = '/Users/jschapMac/Documents/HydrologyData/Livneh/MetNC';

% Directory where clipped forcing files should be saved
forcingsavedir = '/Users/jschapMac/Desktop/Tuolumne/ClippedForcings';

% Name and location of basin shapefile
shpname = '/Users/jschapMac/Desktop/Tuolumne/Shapefiles/upper_tuolumne_wgs.shp';

% Number of forcings in the daily CONUS met. forcing file
numforcings = 4;

% Beginning and ending years of simulation (must be included in the daily CONUS
% met. forcing file)
beginyear = 1985;
endyear = 2011;

% Number of decimal points of precision to use for forcing file names
grid_decimal = 4;

%% Make the mask

tic

addpath(forcingpath)
metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');
    
% Convert lon coords if they use E/W, instead of absolute value system
if max(metlon)>180
    metlon = metlon - 360;
elseif min(metlon)<-180
    % lon = lon some other conversion
end

basin = shaperead(shpname);

if exist('forcmask.mat','file') % only performs masking if necessary
    load forcmask.mat
    disp('Mask loaded.')
else
    mask = NaN(922, 444);
    for i=1:922
        for j=1:444
            mask(i,j) = inpolygon(metlon(i),metlat(j),basin.X,basin.Y);
        end
    end
    save('forcmask.mat', 'mask')
    disp('Mask generated.')
end

[ind1, ind2] = find(mask);

ncells = length(ind1);
nyears = endyear - beginyear + 1;
t_ind = 1;
cum_days = 0;

toc

%% Use the mask to extract forcing variables
disp('Extracting forcing variables')

for t = beginyear:endyear
    
    if t==beginyear, tic, end
    prec = ncread(['prec.' num2str(t) '.nc'], 'prec');
    tmax = ncread(['tmax.' num2str(t) '.nc'], 'tmax');
    tmin = ncread(['tmin.' num2str(t) '.nc'], 'tmin');
    wind = ncread(['wind.' num2str(t) '.nc'], 'wind');
       
    info = ncinfo(['prec.' num2str(t) '.nc']);
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells);
    
    for i=1:ncells
        data(:,1,i) = prec(ind1(i), ind2(i), :);
        data(:,2,i) = tmin(ind1(i), ind2(i), :);
        data(:,3,i) = tmax(ind1(i), ind2(i), :);
        data(:,4,i) = wind(ind1(i), ind2(i), :);
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
for i=1:ncells     
    filename = ['data_' num2str(metlat(ind2(i)),fstring) '_' num2str(metlon(ind1(i)),fstring)];
    dlmwrite(fullfile(forcingsavedir, filename), data_cum(:,:,i))  
end
display(['Forcing data saved to ' forcingsavedir])

%% Extract soils data

% Load the soils data

soils = load(fullfile(soilpath, soilname));

slat = soils(:,3);
slon = soils(:,4);

% Go through each row of soils, check if it is in the study domain.
% If it is, add that row to soilsnew.

disp('Clipping soils data')
ind = 1;
for i=1:length(soils)
    indomain = inpolygon(slon(i),slat(i), basin.X, basin.Y);
    if indomain
        soilsclip(ind,:) = soils(i,:);
        ind = ind + 1;
    end
end

fstring = ['%.' num2str(grid_decimal) 'f'];
fspec = ['%d %.0f ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d'];
dlmwrite(fullfile(soilsavedir, 'soils.TUOLUMNE'), soilsclip, ...
    'precision',fspec,'delimiter','')
display(['Soils data saved to ' soilsavedir])

% Note: the delimiter and the format spec must be specified precisely as
% the Stehekin example from the VIC website in order to avoid the error
% about CELL LATITUDE not being found/for VIC to successfully read the soil
% parameter file.
