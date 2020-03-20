% Subsets forcings to a particular basin
%
% Written 3/9/2020 JRS

function temp = subset_forcings(force_in, force_out, beginyear, endyear, grid_decimal, numforcings, maskname)

mkdir(force_out)
disp(['Created directory ', force_out, ' for outputs']);

[basin_mask, ~, lon, lat] = geotiffread2(maskname);
[mask_lonv, mask_latv, ~] = grid2xyz(lon', lat', basin_mask);
ncells = length(mask_lonv);

metlat = ncread(fullfile(force_in, ['prec.' num2str(beginyear) '.nc']), 'lat');
metlon = ncread(fullfile(force_in, ['prec.' num2str(endyear) '.nc']), 'lon');
metlon = metlon - 360;
        
lat_ind = zeros(ncells,1);
lon_ind = zeros(ncells,1);
for k=1:ncells
        [~, lat_ind(k)] = min(abs(mask_latv(k) - metlat));
        [~, lon_ind(k)] = min(abs(mask_lonv(k) - metlon));
end

%% Extract forcing variables

disp('Extracting forcing variables')

nyears = endyear - beginyear + 1;
t_ind = 1;
cum_days = 0;

for t = beginyear:endyear
    
    if t==beginyear, tic; end
    
    prec = ncread(fullfile(force_in, ['prec.' num2str(t) '.nc']), 'prec');
    tmax = ncread(fullfile(force_in, ['tmax.' num2str(t) '.nc']), 'tmax');
    tmin = ncread(fullfile(force_in, ['tmin.' num2str(t) '.nc']), 'tmin');
    wind = ncread(fullfile(force_in, ['wind.' num2str(t) '.nc']), 'wind');
    
    prec = permute(prec, [2,1,3]);
    tmax = permute(tmax, [2,1,3]);
    tmin = permute(tmin, [2,1,3]);
    wind = permute(wind, [2,1,3]);
       
    % test plot to check the data were loaded properly
%     tmax_map = tmax(:,:,1);
%     figure
%     plotraster([min(metlon) max(metlon)], [min(metlat) max(metlat)], tmax_map, 'tmax', 'Lon','Lat')
    
    info = ncinfo(fullfile(force_in, ['prec.' num2str(t) '.nc']));
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells);

    for k=1:ncells             
        data(:,1,k) = prec(lat_ind(k),lon_ind(k), :);        
        data(:,2,k) = tmin(lat_ind(k),lon_ind(k), :);
        data(:,3,k) = tmax(lat_ind(k),lon_ind(k), :);
        data(:,4,k) = wind(lat_ind(k),lon_ind(k), :);
    end
       
    if t_ind~=1
        data_cum = vertcat(data_cum, data);
    else
        data_cum = data;
    end

    t_ind = t_ind + 1;
    cum_days = size(prec,3) + cum_days;
    
    if t==beginyear 
        disp(['About ' num2str(toc*nyears/60) ' minutes remaining.'])
    end
    disp(['Finished processing year ' num2str(t)])
    
end

%% Save clipped met. forcings

fstring = ['%.' num2str(grid_decimal) 'f'];
for k=1:ncells     
    filename = ['data_' num2str(metlat(lat_ind(k)),fstring) '_' num2str(metlon(lon_ind(k)),fstring)];
    dlmwrite(fullfile(force_out, filename), data_cum(:,:,k), ' ')
end
disp(['Forcing data saved to ' force_out])

temp = data;

return

