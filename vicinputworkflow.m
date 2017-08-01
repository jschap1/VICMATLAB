% Generic workflow for creating VIC input files

soils = load('/Users/jschapMac/Desktop/UpperKaweah/soils.KAWEAH');
slat = soils(:,3);
slon = soils(:,4);

[X Y] = meshgrid(metlat, metlon);

metcoords = NaN(922*444,2);
for i=1:922*444
    metcoords(i,:) = [X(i) Y(i)];
end

i=1;
[a,ind] = ismember([slat(i) slon(i)], metcoords);
% if five decimals of precision are used for the soils data, then the
% coordinates of the Livneh data match them.

% Extract the forcing field values for the cells whose coordinates match
% those of the soils data.



% 36.7188 -119.3438

indlat = metcoords(:,1)>36.7 & metcoords(:,1)<36.72;
indlatlon = (metcoords(:,1)>36.7 & metcoords(:,1)<36.72) & metcoords(:,2)>-119.4 & metcoords(:,2)<-119.3;
metcoords(indlatlon)

% soils result
% 36.7188000000000  -119.3438000000000
% 36.7187500000000  -119.3437500000000

% metcoords result:
% 36.7187500000000  -119.3437500000000

%% Specify inputs

shpname = '/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.shp';
savedir = '/Users/jschapMac/Desktop/UpperKaweah/ClippedForcings2';
numforcings = 4;
beginyear = 2000;
endyear = 2007;
forcingdir = '/Users/jschapMac/Desktop/UpperKaweah/MetNC';

%% Make mask

% Makes a mask of the met. forcing dataset over the CONUS. The extent of
% the mask is the basin of interest.

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

basin = shaperead(shpname);

if exist('mask.mat','file') % only performs masking if necessary
    load mask.mat
    disp('Mask loaded.')
else
    mask = NaN(922, 444);
    for i=1:922
        for j=1:444
            mask(i,j) = inpolygon(metlon(i),metlat(j),basin.X,basin.Y);
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

extract_forcing_over_region
