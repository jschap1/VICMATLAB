% This is better as a script

cd '/Users/jschapMac/Desktop/UpperKaweah'; % directory containing basin shapefile
kaweah = shaperead('upper_kaweah.shp'); % name of basin shapefile
numforcings = 4;

% ncinfo('MetNC/prec.2007.nc')
lat = ncread('MetNC/prec.2007.nc', 'lat'); % load forcing netCDF files
lon = ncread('MetNC/prec.2007.nc', 'lon');
time = ncread('MetNC/prec.2007.nc', 'time');
prec = ncread('MetNC/prec.2007.nc', 'prec');
tmax = ncread('MetNC/tmax.2007.nc', 'tmax');
tmin = ncread('MetNC/tmin.2007.nc', 'tmin');
wind = ncread('MetNC/wind.2007.nc', 'wind');

% Convert lon coords if they use E/W, instead of absolute value system
if max(lon)>180
    lon = lon - 360;
elseif min(lon)<-180
    % lon = lon some other conversion
end

if exist('mask.mat','file') % only performs masking if necessary
    load mask.mat
else
    mask = NaN(922, 444);
    for i=1:922
        for j=1:444
            mask(i,j) = inpolygon(lon(i),lat(j),kaweah.X,kaweah.Y);
        end
    end
    save('mask.mat', 'mask')
end

[ind1, ind2] = find(mask);

ncells = length(ind1);
for i=1:ncells
    data = NaN(365,numforcings);
    data(:,1) = prec(ind1(i), ind2(i), :);
    data(:,2) = tmin(ind1(i), ind2(i), :);
    data(:,3) = tmax(ind1(i), ind2(i), :);
    data(:,4) = wind(ind1(i), ind2(i), :);
    filename = ['data_' num2str(lat(ind2(i))) '_' num2str(lon(ind1(i)))];
    dlmwrite(fullfile('ClippedForcings', [filename '.txt']), data)
end

% figure, hold on 
% plot(kaweah.X,kaweah.Y,'LineWidth',2)
% plot(lon(ind1), lat(ind2), '*')
% hold off