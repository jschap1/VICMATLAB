% Clips soil parameter file to basin extent
%
% Adapted from vicinputworkflow 3/12/2019 JRS
%
% INPUTS
% Soil parameter file encompassing a region at least as large as the basin
% List of coordinates in basin
% 
% OUTPUTS
% Clipped soil parameter file
%
% Double check that this is working properly, especially for elevation
% 4/4/2019

cd '/Volumes/HD3/SWOTDA/Data/IRB/'

soilpath = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_0';
soilname = 'soils_3L_MERIT.txt'; 
soilsavedir = './VIC';
grid_decimal = 5; % precision used in forcing filenames

% Use r.out.xyz to generate this from the basin mask raster
maskxyz = dlmread('basincoords.txt', '|');
masklon = maskxyz(maskxyz(:,3) == 1,1);
masklat = maskxyz(maskxyz(:,3) == 1,2);
ncells = length(masklon);

soils = load(fullfile(soilpath, soilname));
slat = soils(:,3);
slon = soils(:,4);

soils(:,1) = 0; % turn off all cells in the soil parameter file

% turn on the grid cells in the soils file that match the grid cells in the
% mask

resolution = 1/16;

for k=1:ncells
    
    nearby_lats = find(abs(slat-masklat(k)) < 0.5*resolution);
    nearby_lons = find(abs(slon-masklon(k)) < 0.5*resolution);
    
    ind = find(ismember(nearby_lats, nearby_lons));
    cellnumber = nearby_lats(ind);
%     slat(cellnumber)
%     slon(cellnumber)
    soils(cellnumber,1) = 1;
    
    % report progres
    if mod(k, 1e3) == 0
        disp(k)
    end
    
end

soils_copy = soils;
soils(soils(:,1) == 0,:) = []; % include this line to reduce file size

if size(soils, 1) < ncells
    disp('Some mask grid cells were not in the soil parameter file')
    disp([num2str(ncells - size(soils, 1)) ' to be exact'])
end

% Write out
fstring = ['%.' num2str(grid_decimal) 'f'];

% For the Livneh soil parameter file
% fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];

% For the 2-layer soil parameter file from HWSD
% fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %d %d %.3f %.3f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d\n'];

% This is used VIC-3L w no optional variables (54 columns in the soil parameter file).
fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %.2f\n'];

fID = fopen(fullfile(soilsavedir, 'soils.SB'),'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved to ' soilsavedir])