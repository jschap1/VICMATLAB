% Convert forcing
%
% Converts ASCII forcing data inputs to NetCDF forcings
% Written 4/14/2020 JRS
% Based on code written by Dongyue Li

function convert_forcing(indir, prefix, outname, precision, start_date, end_date, nt_per_day)

disp('reading forcing file names')
fnames = dir(fullfile(indir, [prefix '*']));

% for testing
% fnames = fnames(1:3);

% info = ncinfo(fullfile(indir, fnames(1).name));
disp('read forcing file names')

ncells = length(fnames);
tmpforc = dlmread(fullfile(indir, fnames(1).name)); 
nsteps = size(tmpforc,1);
nvars = size(tmpforc,2);

% Adapted from LoadVICResults
gridcells = cell(ncells, 1);
for k=1:ncells
    tmpstring = fnames(k).name;
    tmpstring = strrep(tmpstring,'.','_');
    gridcells{k} = tmpstring;
end

% precision = 5;
[lat, lon] = GetCoords(gridcells, precision, prefix);
% resolution = max(abs(lat(2)-lat(1)), abs(lon(2)-lon(1)));

% lon = min(lon):resolution:max(lon);
% lat = min(lat):resolution:max(lat);
% nlat = length(lat);
% nlon = length(lon);

sample_forc = dlmread(fullfile(indir, fnames(1).name));
nt = size(sample_forc,1);
% nt_per_day = 24; % hard-coded (assumes hourly forcing data)
% ndays = nt/nt_per_day;
% start_date = datetime(1980, 1, 1, 0, 0, 0);
% end_date = datetime(2011, 12, 31, 23, 0, 0);
timevector = start_date:hours(1):end_date;
timevector = timevector';

mask=ones(length(lon),length(lat));
mask = single(mask);

nlon = length(unique(lon));
nlat = length(unique(lat));

% Prepare inputs for write_netcdf_forcing
info.lon = lon;
info.lat = lat;
info.mask = mask;

%% Checks

if length(timevector) ~= nt
    error('check start and end dates')
end

if month(start_date) ~= 1 || day(start_date) ~= 1
    error('start date must be Jan. 1')
end

if month(end_date) ~= 12 || day(end_date) ~= 31
    error('end date must be Dec. 31')
end

%% Main loop

disp('beginning main loop')

yy = year(start_date);
ndays_in_year = yeardays(yy);
t1 = 1;
t2 = ndays_in_year*nt_per_day;
    
while yy <= year(end_date)
    
    ndays_in_year = yeardays(yy);
    
    disp(['current year is ' num2str(yy)]);
    
    temperature = NaN(ndays_in_year*nt_per_day, ncells);
    precipitation = NaN(ndays_in_year*nt_per_day, ncells);
    pressure = NaN(ndays_in_year*nt_per_day, ncells);
    shortwave = NaN(ndays_in_year*nt_per_day, ncells);
    longwave = NaN(ndays_in_year*nt_per_day, ncells);
    vp = NaN(ndays_in_year*nt_per_day, ncells);
    wind = NaN(ndays_in_year*nt_per_day, ncells);
    
    for k=1:ncells
        
        forc = dlmread(fullfile(indir, fnames(k).name));

        % indexes are based on the order of the variables in the ASCII
        % forcing file
        temperature(:, k) = forc(t1:t2,2);
        precipitation(:, k) = forc(t1:t2,1);
        pressure(:, k) = forc(t1:t2,6);
        shortwave(:, k) = forc(t1:t2,3);
        longwave(:, k) = forc(t1:t2,4);
        vp(:, k) = forc(t1:t2,7);
        wind(:, k) = forc(t1:t2,8);

        disp(k)
    end
    
    temperature_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    precipitation_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    pressure_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    shortwave_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    longwave_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    vp_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    wind_map = NaN(nlon, nlat, ndays_in_year*nt_per_day);
    for t=1:ndays_in_year*nt_per_day
        temperature_map(:,:,t) = fliplr(xyz2grid(lon, lat, temperature(t,:)')');
        precipitation_map(:,:,t) = fliplr(xyz2grid(lon, lat, precipitation(t,:)')');
        pressure_map(:,:,t) = fliplr(xyz2grid(lon, lat, pressure(t,:)')');
        shortwave_map(:,:,t) = fliplr(xyz2grid(lon, lat, shortwave(t,:)')');
        longwave_map(:,:,t) = fliplr(xyz2grid(lon, lat, longwave(t,:)')');
        vp_map(:,:,t) = fliplr(xyz2grid(lon, lat, vp(t,:)')');
        wind_map(:,:,t) = fliplr(xyz2grid(lon, lat, wind(t,:)')');
    end
        
%     figure, plotraster(lon, lat, temperature_map(:,:,1), 'temp')
    
    info.ndays = ndays_in_year;
    info.nt_per_day = nt_per_day;
    info.outname = [outname '_' num2str(yy) '.nc'];
    write_netcdf_forcing(temperature_map, precipitation_map, pressure_map, ...
        shortwave_map, longwave_map, vp_map, wind_map, info)

    yy = yy + 1; % iterate/move on to next year
    ndays_in_year = yeardays(yy);
    t1 = t2 + 1;
    t2 = t1 + ndays_in_year*nt_per_day - 1;
    
end
