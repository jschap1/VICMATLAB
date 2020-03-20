%% Specify inputs

shpname = '/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.shp';
savedir = '/Users/jschapMac/Desktop/UpperKaweah/ClippedForcings2';
numforcings = 4;
beginyear = 2000;
endyear = 2007;
forcingdir = '/Users/jschapMac/Desktop/UpperKaweah/MetNC';

%% Read in Livneh data coordinates
tic

addpath(forcingdir)
metlat = ncread(['prec.' num2str(beginyear) '.nc'], 'lat');
metlon = ncread(['prec.' num2str(beginyear) '.nc'], 'lon');
    
% Convert lon coords if they use E/W, instead of absolute value system
if max(metlon)>180
    metlon = metlon - 360;
elseif min(metlon)<-180
    % lon = lon some other conversion
end

ncells = 103;
nyears = endyear - beginyear + 1;
t_ind = 1;
cum_days = 0;

toc

%% Use the soil parameter file to extract forcing variables
disp('Extracting forcing variables')

soils = load('/Users/jschapMac/Desktop/UpperKaweah/soils.KAWEAH');
slat = soils(:,3);
slon = soils(:,4);

% Get the index of the Livneh grid that matches the soil parameter file lat
% lon values

[X Y] = meshgrid(metlat, metlon);

metcoords = NaN(922*444,2);
for i=1:922*444
    metcoords(i,:) = [X(i) Y(i)];
end

indlat = NaN(103,1);
indlon = NaN(103,1);
for i=1:103
    [a,indlat(i)] = ismember(slat(i), metcoords(:,1));
    [a,indlon(i)] = ismember(slon(i), metcoords(:,2));
end



% ind = NaN(103,2);
% for i=1:103
%     [a,ind(i,:)] = ismember([slat(i) slon(i)], metcoords);
% end

prec(indlon(i), indlat(i), :)

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
    filename = ['data_' num2str(metlat(ind2(i))) '_' num2str(metlon(ind1(i)))];
    dlmwrite(fullfile(savedir, [filename '.txt']), data_cum(:,:,i))  
end