% Given a soil parameter file for the entire CONUS, selects the rows
% corresponding to cells in the model domain, and creates a subsetted soil
% parameter file containing only data for the grid cells in the study
% domain.

%% Specify inputs

maskname = '/Users/jschapMac/Desktop/UpperKaweah/mask.mat';
savedir = '/Users/jschapMac/Desktop/UpperKaweah/ClippedSoils';
savename = 'soils.KAWEAH';

% name of full soil parameter file
soilpath = '/Users/jschapMac/Documents/HydrologyData/VICParametersCONUS';
soilname = 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt'; 

% loads lat and lon for study domain
shpname = '/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.shp';
basin = shaperead(shpname);
domainlat = basin.Y;
domainlon = basin.X;

%% Load NLDAS soil parameter file

soils = load(fullfile(soilpath, soilname));

lat = soils(:,3);
lon = soils(:,4);

%% Make a soil mask

% Go through each row of soils, check if it is in the study domain.
% If it is, add that row to soilsnew.

ind = 1;
for i=1:length(soils)
    indomain = inpolygon(lon(i),lat(i), domainlon, domainlat);
    if indomain
        soilsnew(ind,:) = soils(i,:);
        ind = ind + 1;
    end
end

dlmwrite(fullfile(savedir, savename), soilsnew)
