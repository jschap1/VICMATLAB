% Subsets forcings to a particular basin
%
% Written 3/9/2020 JRS
%
% Supercedes subset_forcings
% Code has been modified to use MUCH less RAM
% 
% Default values for optional arguments
% numforcings = 4;
% grid_decimal = 5;
% has_lots_of_ram = 0;
%
% Requires PREC, TMIN, TMAX, WIND forcings in NetCDF files, one per year

function data_cum = subset_forcings(indir, outdir, beginyear, endyear, maskname, varargin)

% Set default values for optional arguments

numvarargs = length(varargin);
if numvarargs > 3
    error('subset_forcings requires at most three optional inputs');
end

optargs = {4, 5, 0};
optargs(1:numvarargs) = varargin;

[numforcings, grid_decimal, has_lots_of_ram] = optargs{:};

mkdir(outdir)
disp(['Created directory ', outdir, ' for outputs']);

% Load basin mask
[basin_mask, ~, lon, lat] = geotiffread2(maskname);

if isa(basin_mask, 'single')
    disp('basin mask is single')
    disp('converting zeros to nans for subsetting')
    basin_mask = double(basin_mask);
    basin_mask(basin_mask==0) = NaN;
end

[mask_lonv, mask_latv, ~] = grid2xyz(lon', lat', basin_mask);
ncells = length(mask_lonv);

% Identify grid cells to keep
metlat = ncread(fullfile(indir, ['prec.' num2str(beginyear) '.nc']), 'lat');
metlon = ncread(fullfile(indir, ['prec.' num2str(endyear) '.nc']), 'lon');
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
    
    prec = ncread(fullfile(indir, ['prec.' num2str(t) '.nc']), 'prec');
    tmax = ncread(fullfile(indir, ['tmax.' num2str(t) '.nc']), 'tmax');
    tmin = ncread(fullfile(indir, ['tmin.' num2str(t) '.nc']), 'tmin');
    wind = ncread(fullfile(indir, ['wind.' num2str(t) '.nc']), 'wind');
    % each of these arrays requires about 1.2 GB of storage
    
    prec = permute(prec, [2,1,3]);
    tmax = permute(tmax, [2,1,3]);
    tmin = permute(tmin, [2,1,3]);
    wind = permute(wind, [2,1,3]);
       
    % test plot to check the data were loaded properly
%     tmax_map = tmax(:,:,1);
%     figure
%     plotraster([min(metlon) max(metlon)], [min(metlat) max(metlat)], tmax_map, 'tmax', 'Lon','Lat')
    
    info = ncinfo(fullfile(indir, ['prec.' num2str(t) '.nc']));
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells);

    for k=1:ncells             
        data(:,1,k) = prec(lat_ind(k),lon_ind(k), :);        
        data(:,2,k) = tmin(lat_ind(k),lon_ind(k), :);
        data(:,3,k) = tmax(lat_ind(k),lon_ind(k), :);
        data(:,4,k) = wind(lat_ind(k),lon_ind(k), :);
    end
       
    % data_cum grows on each loop iteration
    % unless the computer has a ton of RAM, this is going to get out of
    % hand, so there is an option to save the data to disk and read it back
    % in, in a separate step
    
    if has_lots_of_ram
        cum_days = size(prec,3) + cum_days;
        if t_ind~=1 
            data_cum = vertcat(data_cum, data);
        else
            data_cum = data;
        end
    else
        disp('save this year''s processed data to file')
        
        savename = fullfile(outdir, ['clipped_forcing_data_year_' num2str(t) '.mat']);
        save(savename, 'data');
        
        cum_days = size(prec,3) + cum_days;
        clearvars('prec','tmax','tmin','wind');
    end
    
    t_ind = t_ind + 1;

    if t==beginyear 
        disp(['About ' num2str(toc*nyears/60) ' minutes remaining.'])
    end
    disp(['Finished processing year ' num2str(t)])
    
end


%% Combine forcings together, if necessary
if ~has_lots_of_ram
    t_ind = 1;
    for t=beginyear:endyear
        savename = fullfile(outdir, ['clipped_forcing_data_year_' num2str(t) '.mat']);
        load(savename, 'data');
        if t_ind==1
            data_cum = data;
        else
            data_cum = vertcat(data_cum, data);
        end
        t_ind = t_ind + 1;
    end
end

%% Save clipped met. forcings
save_forcings(metlat, lat_ind, metlon, lon_ind, outdir, data_cum, grid_decimal);

end


% fstring = ['%.' num2str(grid_decimal) 'f'];
% for k=1:ncells     
%     filename = ['data_' num2str(metlat(lat_ind(k)),fstring) '_' num2str(metlon(lon_ind(k)),fstring)];
%     dlmwrite(fullfile(force_out, filename), data_cum(:,:,k), ' ')
% end
% disp(['Forcing data saved to ' force_out])