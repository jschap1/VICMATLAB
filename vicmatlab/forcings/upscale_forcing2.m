% Rescales forcings to a different resolution
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

function data_cum = upscale_forcing2(indir, outdir, beginyear, endyear, newres, oldres, varargin)

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

% Identify grid cells
metlat = ncread(fullfile(indir, ['prec.' num2str(beginyear) '.nc']), 'lat');
metlon = ncread(fullfile(indir, ['prec.' num2str(endyear) '.nc']), 'lon');
metlon = metlon - 360;

% f = newres/oldres; % scaling factor

% m1 = length(metlat);
% n1 = length(metlon);

% new dimensions
% m2 = m1/f;
% n2 = floor(n1/f);

% Upscale (or downscale) (method from Walter Roberson in MATLAB Answers)
% minlat = min(metlat);
% maxlat = max(metlat);
% minlon = min(metlon);
% maxlon = max(metlon);
% oldlats = linspace(minlat, maxlat, m1);
% oldlons = linspace(minlon, maxlon, n1)';
% newlats = linspace(minlat, maxlat, m2);
% newlons = linspace(minlon, maxlon, n2)';

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
       
    if t==beginyear % resample maskfile to get output grid
        
        % test plot to check the data were loaded properly
        tmax_map = tmax(:,:,1);
        figure
        plotraster([min(metlon) max(metlon)], [min(metlat) max(metlat)], tmax_map, 'tmax', 'Lon','Lat')
        
        forcmask = tmax_map;
        forcmask(~isnan(forcmask)) = 1;
        figure
        plotraster([min(metlon) max(metlon)], [min(metlat) max(metlat)], forcmask, 'mask', 'Lon','Lat')
        
%         forcmask2 = interp2(oldlats, oldlons, forcmask', newlats, newlons, 'linear');
        [forcmask2, newlons, newlats] = upscale_raster(forcmask, metlon, metlat, newres, oldres, 'linear');
        
%         forcmask2 = forcmask2';
        forcmask2(isnan(forcmask2))=0;
        forcmask3 = logical(forcmask2(:));
        
        ncells1 = nansum(forcmask(:));
        ncells2 = nansum(forcmask2(:));
        
    end
    
    info = ncinfo(fullfile(indir, ['prec.' num2str(t) '.nc']));
    ndays = info.Dimensions(1).Length; % get number of days in the year
    data = NaN(ndays, numforcings, ncells2);

    [m2, n2] = size(forcmask2);
    
    prec2 = interp3(metlat', metlon, 1:ndays, permute(prec, [2,1,3]), newlats, newlons, 1:ndays, 'linear');
    prec2 = permute(prec2, [2,1,3]);
    prec2(isnan(prec2)) = 0;
    prec3 = reshape(prec2, n2*m2, ndays);
    prec3 = prec3(logical(forcmask3),:);
    data(:,1,:) = permute(prec3, [2,1]);
    
    tmax2 = interp3(metlat', metlon, 1:ndays, permute(tmax, [2,1,3]), newlats, newlons, 1:ndays, 'linear');
    tmax2 = permute(tmax2, [2,1,3]);
    tmax2(isnan(tmax2)) = 0;
    tmax3 = reshape(tmax2, n2*m2, ndays);
    tmax3 = tmax3(logical(forcmask3),:);
    data(:,2,:) = permute(tmax3, [2,1]);    
    
    tmin2 = interp3(metlat', metlon, 1:ndays, permute(tmin, [2,1,3]), newlats, newlons, 1:ndays, 'linear');
    tmin2 = permute(tmin2, [2,1,3]);
    tmin2(isnan(tmin2)) = 0;
    tmin3 = reshape(tmin2, n2*m2, ndays);
    tmin3 = tmin3(logical(forcmask3),:);
    data(:,3,:) = permute(tmin3, [2,1]);

    wind2 = interp3(metlat', metlon, 1:ndays, permute(wind, [2,1,3]), newlats, newlons, 1:ndays, 'linear');
    wind2 = permute(wind2, [2,1,3]);
    wind2(isnan(wind2)) = 0;
    wind3 = reshape(wind2, n2*m2, ndays);
    wind3 = wind3(logical(forcmask3),:);
    data(:,4,:) = permute(wind3, [2,1]);
    
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

%% Save upscaled met. forcings

[X,Y] = ndgrid(newlons, newlats);

% [xx1,yy,~] = grid2xyz(X,Y,forcmask2);
[xx,yy,~] = grid2xyz(newlons, newlats,forcmask2);
xx = xx(forcmask3);
yy = yy(forcmask3);

fstring = ['%.' num2str(grid_decimal) 'f'];
for k=1:ncells2     
    filename = ['data_' num2str(yy(k),fstring) '_' num2str(xx(k),fstring)];
    dlmwrite(fullfile(outdir, filename), data_cum(:,:,k), ' ')
end
disp(['Forcing data saved to ' outdir])

end
