% Convert forcing
%
% Converts ASCII forcing data inputs to NetCDF forcings
% Written 4/14/2020 JRS
% Based on code written by Dongyue Li

function convert_forcing(indir, prefix, outname, precision, start_date, end_date, nt_per_day)

delt = hours(24/nt_per_day);
if mod(24/nt_per_day,1) ~=0
    error('nt_per_day must divide 24 evenly')
end

disp('reading forcing file names')
fnames = dir(fullfile(indir, [prefix '*']));

% for testing
% fnames = fnames(1:3);

% info = ncinfo(fullfile(indir, fnames(1).name));
disp('read forcing file names')

ncells = length(fnames);
% tmpforc = dlmread(fullfile(indir, fnames(1).name));
% clearvars tmpforc
% nsteps = size(tmpforc,1);
% nvars = size(tmpforc,2);

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
clearvars sample_forc

% nt_per_day = 24; % hard-coded (assumes hourly forcing data)
% ndays = nt/nt_per_day;
% start_date = datetime(1980, 1, 1, 0, 0, 0);
% end_date = datetime(2011, 12, 31, 23, 0, 0);

timevector = start_date:delt:end_date;
timevector = timevector';

nlon = length(unique(lon));
nlat = length(unique(lat));

% Prepare inputs for write_netcdf_forcing
info.lon = lon;
info.lat = lat;
info.start_year = year(start_date);

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

    
%% Calculate indexes before entering the main loop

yy = year(start_date);
ndays_in_year = yeardays(yy);
yearvect = year(start_date):year(end_date);
nyears = length(yearvect);
t1 = zeros(nyears,1);
t2 = zeros(nyears,1);
t1(1) = 1;
t2(1) = ndays_in_year*nt_per_day;

for k=2:nyears
    ndays_in_year = yeardays(yearvect(k));
    t1(k) = t2(k-1) + 1;
    t2(k) = t1(k) + ndays_in_year*nt_per_day - 1;    
end

%% Main loop

disp('beginning main loop')

% while yy <= year(end_date)

% c = parcluster('local');
% nw = c.NumWorkers;
% parpool(nw-1)
% parfor j=1:nyears
for j=1:nyears    
    
    yy = yearvect(j);
    ndays_in_year = yeardays(yy);
    
    disp(['current year is ' num2str(yy)]);
    
    temperature = NaN(ndays_in_year*nt_per_day, ncells);
    precipitation = NaN(ndays_in_year*nt_per_day, ncells);
    pressure = NaN(ndays_in_year*nt_per_day, ncells);
    shortwave = NaN(ndays_in_year*nt_per_day, ncells);
    longwave = NaN(ndays_in_year*nt_per_day, ncells);
    vp = NaN(ndays_in_year*nt_per_day, ncells);
    wind = NaN(ndays_in_year*nt_per_day, ncells);
    
    tic
    for k=1:ncells
        
        forc = dlmread(fullfile(indir, fnames(k).name));

        % indexes are based on the order of the variables in the ASCII
        % forcing file
        temperature(:, k) = forc(t1(j):t2(j),2);
        precipitation(:, k) = forc(t1(j):t2(j),1);
        pressure(:, k) = forc(t1(j):t2(j),6);
        shortwave(:, k) = forc(t1(j):t2(j),3);
        longwave(:, k) = forc(t1(j):t2(j),4);
        vp(:, k) = forc(t1(j):t2(j),7);
        wind(:, k) = forc(t1(j):t2(j),8);

        if mod(k,100) == 0
            disp(k)
        end
        
        if k==10
            toc
            disp(['Estimated time remaining: ' num2str((toc*ncells/10/3600)), ' hours'])
        end
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
    
    clearvars temperature precipitation pressure shortwave longwave vp wind

    info.ndays = ndays_in_year;
    info.nt_per_day = nt_per_day;
    info.outname = [outname '.' num2str(yy) '.nc'];
    info.year = yy;
    
%     save('/hdd/ESSD/data/stehekin/write_nc_forc_inputs.mat','temperature_map', 'precipitation_map', 'pressure_map', ...
%         'shortwave_map', 'longwave_map', 'vp_map', 'wind_map', 'info')
    
    write_netcdf_forcing(temperature_map, precipitation_map, pressure_map, ...
        shortwave_map, longwave_map, vp_map, wind_map, info)

%     yy = yy + 1; % iterate/move on to next year
%     ndays_in_year = yeardays(yy);
%     t1 = t2 + 1;
%     t2 = t1 + ndays_in_year*nt_per_day;
    % if I move t1, t2, ndays_in_year out of the loop, then I can make it a
    % parfor loop.
    %
    % however, there is a RAM limitation. Will most likely need to use
    % 8-bit integers to store information to save space
    
end

end