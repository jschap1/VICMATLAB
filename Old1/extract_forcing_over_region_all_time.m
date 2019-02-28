% Finds the precise forcing data that applies to the region of interest,
% from a larger forcing dataset. Must run from directory containing your
% basin shapefiles and forcing data.
%
% It may take a long time to run (hours to days) depending on the number of
% years and the size of the modeling domain.

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
for i=1:ncells
    
    t_ind = 1;
    cum_days = 0; % needed to keep track of the length of "data"
    
    for t = beginyear:endyear
        
        % Amend this to just read a small part of the nc file, in
        % particular the part containing the basin of interest.
        
        if t==beginyear && i==1, tic, end
        prec = ncread(['MetNC/prec.' num2str(t) '.nc'], 'prec');
        tmax = ncread(['MetNC/tmax.' num2str(t) '.nc'], 'tmax');
        tmin = ncread(['MetNC/tmin.' num2str(t) '.nc'], 'tmin');
        wind = ncread(['MetNC/wind.' num2str(t) '.nc'], 'wind');
        
        % Note: in order to pre-allocate "data", would need to know the
        % number of days in advance, instead of reading from prec.
        data(:,1,t_ind) = prec(ind1(i), ind2(i), :);
        data(:,2,t_ind) = tmin(ind1(i), ind2(i), :);
        data(:,3,t_ind) = tmax(ind1(i), ind2(i), :);
        data(:,4,t_ind) = wind(ind1(i), ind2(i), :);
        
        t_ind = t_ind + 1;
        
        cum_days = size(prec,3) + cum_days;
        
        if t==beginyear && i==1, 
            disp(['About ' num2str(toc*ncells*nyears/60) ...
                ' minutes remaining.'])
        end
        
    end
    
    savedata = reshape(data,cum_days, numforcings);
        
    filename = ['data_' num2str(lat(ind2(i))) '_' num2str(lon(ind1(i)))];
    dlmwrite(fullfile('ClippedForcings', [filename '.txt']), savedata)
    
end