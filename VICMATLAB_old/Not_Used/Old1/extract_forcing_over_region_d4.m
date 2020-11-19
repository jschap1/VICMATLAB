% Testing to see if I can load only a portion of the nc files

tic
cd '/Users/jschapMac/Desktop/UpperKaweah'; % directory containing basin shapefile
kaweah = shaperead('upper_kaweah.shp'); % name of basin shapefile
numforcings = 4;
beginyear = 2006;
endyear = 2007;

lat = ncread(['MetNC/prec.' num2str(beginyear) '.nc'], 'lat');
lon = ncread(['MetNC/prec.' num2str(beginyear) '.nc'], 'lon');
    
% Convert lon coords if they use E/W, instead of absolute value system
if max(lon)>180
    lon = lon - 360;
elseif min(lon)<-180
    % lon = lon some other conversion
end

if exist('mask.mat','file') % only performs masking if necessary
    load mask.mat
    disp('Mask loaded.')
else
    mask = NaN(922, 444);
    for i=1:922
        for j=1:444
            mask(i,j) = inpolygon(lon(i),lat(j),kaweah.X,kaweah.Y);
        end
    end
    save('mask.mat', 'mask')
    disp('Mask generated.')
end

[ind1, ind2] = find(mask);

ncells = length(ind1);
nyears = endyear - beginyear + 1;
toc

disp('Extracting forcing variables')

for t = beginyear:endyear

    info = ncinfo(['MetNC/prec.' num2str(t) '.nc']);
    ndays = info.Dimensions(1).Length;
    start = [min(ind1),min(ind2),1];
    count = [max(ind1),max(ind2),ndays]; % lon, lat, time

%     prec = ncread(['MetNC/prec.' num2str(t) '.nc'], 'prec');
%     tmax = ncread(['MetNC/tmax.' num2str(t) '.nc'], 'tmax');
%     tmin = ncread(['MetNC/tmin.' num2str(t) '.nc'], 'tmin');
%     wind = ncread(['MetNC/wind.' num2str(t) '.nc'], 'wind');
    
    prec = ncread(['MetNC/prec.' num2str(t) '.nc'], 'prec', start, count);
    tmax = ncread(['MetNC/tmax.' num2str(t) '.nc'], 'tmax', start, count);
    tmin = ncread(['MetNC/tmin.' num2str(t) '.nc'], 'tmin', start, count);
    wind = ncread(['MetNC/wind.' num2str(t) '.nc'], 'wind', start, count);
    
    for i=1:ncells
        data(:,1,i,t_ind) = prec(ind1(i), ind2(i), :);
        data(:,2,i,t_ind) = tmin(ind1(i), ind2(i), :);
        data(:,3,i,t_ind) = tmax(ind1(i), ind2(i), :);
        data(:,4,i,t_ind) = wind(ind1(i), ind2(i), :);
    end

    t_ind = t_ind + 1;
    
end

% met. forcing variable statistics for a particular day and grid cell
summary(table(data(:,:,1,1)));

%%
% Check whether the numbers are correct.

% data(15,4,17,2) should be the wind speed on day 15 for pixel 17 in year
% 2007.

wind_check = ncread(['MetNC/wind.' num2str(2007) '.nc'], 'wind');

% pixel 17 has the following lon/lat:
lon(ind1(17))
lat(ind2(17))
% (-119.4062, 36.3438)

u_check = wind_check(ind1(17), ind2(17), 15);
u = data(15,4,17,2);

% They are equal, so that is good.
%%
% Check, using the partial method
% Now, they do not match.
% For now, use the full method.

%%


%
%
%
%
%
%
    % Note: in order to pre-allocate "data", would need to know the
    % number of days in advance, instead of reading from prec.
    data(:,1,t_ind) = prec(ind1(i), ind2(i), :);
    data(:,2,t_ind) = tmin(ind1(i), ind2(i), :);
    data(:,3,t_ind) = tmax(ind1(i), ind2(i), :);
    data(:,4,t_ind) = wind(ind1(i), ind2(i), :);

    t_ind = t_ind + 1;

    cum_days = size(prec,3) + cum_days;

end




for i=1:ncells
    
    t_ind = 1;
    cum_days = 0; % needed to keep track of the length of "data"
    clear data
    

    
    savedata = reshape(data,cum_days, numforcings);
        
    filename = ['data_' num2str(lat(ind2(i))) '_' num2str(lon(ind1(i)))];
    savedir = 'ClippedForcingsV2';
    dlmwrite(fullfile(savedir, [filename '.txt']), savedata)
    
end
