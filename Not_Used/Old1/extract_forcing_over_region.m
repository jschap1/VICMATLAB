% Finds the precise forcing data that applies to the region of interest,
% from a larger forcing dataset. 
%
% Must run from directory containing your basin shapefiles and forcing data. 
%
% Takes netCDF (daily Livneh data) as input, and outputs ASCII files 
% compatible with the VIC4 model.
%
% Currently loads the whole NC file, which slows it down unnecessarily
%
% INPUTS
% watershed/basin shapefile
% savedirectory
% nforcings
% beginyear
% endyear
% daily forcing data that includes data over the region of interest in
% NetCDF format. For example, Livneh forcing data.
%
% OUTPUTS
% watershed mask
% clipped forcing files of size (ndays, nforcings) for each grid cell

%% Specify inputs

shpname = '/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.shp';
savedir = '/Users/jschapMac/Desktop/UpperKaweah/ClippedForcings2';
numforcings = 4;
beginyear = 2000;
endyear = 2007;

%% Make the mask

tic

lat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
lon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');
    
% Convert lon coords if they use E/W, instead of absolute value system
if max(lon)>180
    lon = lon - 360;
elseif min(lon)<-180
    % lon = lon some other conversion
end

basin = shaperead(shpname);

if exist('mask.mat','file') % only performs masking if necessary
    load mask.mat
    disp('Mask loaded.')
else
    mask = NaN(922, 444);
    for i=1:922
        for j=1:444
            mask(i,j) = inpolygon(lon(i),lat(j),basin.X,basin.Y);
        end
    end
    save('mask.mat', 'mask')
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

for i=1:ncells     
    filename = ['data_' num2str(lat(ind2(i))) '_' num2str(lon(ind1(i)))];
    dlmwrite(fullfile(savedir, [filename '.txt']), data_cum(:,:,i))  
end