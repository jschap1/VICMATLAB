% For each met. forcing file, finds the closest grid cell from the soil
% parameter file and renames the forcing file to match the grid cell
% coordinates from the soil parameter file
%
% Sometimes it fails to find a nearby forcing file.
%
% Need to be on the correct grid (the same grid as the soil parameter file)
% for the VIC model to run
%
% These lines change the file names by choosing the nearest soil parameter
% grid cell
%
% Also, we need remove the .txt extension in bash so VIC can read the files
% paste -> for file in fluxes*.txt; do mv "$file" "${file%.txt}"; done
%
% This is VERY slow if you run it in the directory with all the forcing
% files. Use relative paths instead. See https://www.mathworks.com/matlabcentral/answers/112086-matlab-slow-when-too-many-files-in-directory

% soils = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_2/soils_3L_MERIT.txt');
% soils = load('./Data/IRB/VIC/MiniDomain2/soils.SB');
% forcingdir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain2/2018-2019_forc_lasttest';
% target_res = 1/16;
% precision = '%3.5f';

function find_closest_soil_gridcell(soils, forcingdir, copyloc, target_res, precision)

% Location to write outputs
if ~exist(copyloc, 'dir')
   mkdir(copyloc)
   disp(['Created directory for outputs: ' copyloc])
else
    disp(['Outputs will be written to: ' copyloc])
end

forcnames = dir(fullfile(forcingdir, 'Forcings_*'));
% soils = load(soilparamfile);

slat = soils(:,3);
slon = soils(:,4);

% cd /Volumes/HD3/SWOTDA/Data/IRB/VIC
% forcnames = dir('./Forc_halfyear/Forcings_*');
% soils = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/soils.SB');
% 
% slat = soils(:,3);
% slon = soils(:,4);

% find closest
% target_res = 1/16;
% precision = '%3.5f';

nmissing = 0;
ncells = length(forcnames);

for k=1:ncells
    
    tmp = strsplit(forcnames(k).name, '_');
    flat = str2double(tmp{2});
    tmp = strsplit(tmp{3}, '.txt');
    flon = str2double(tmp{1});

    nearby_lats = find(abs(slat-flat) < 0.5*target_res);
    nearby_lons = find(abs(slon-flon) < 0.5*target_res);
    ind = find(ismember(nearby_lats, nearby_lons));
    
    if isempty(ind)  
        nmissing = nmissing + 1;     
        continue
    end
    
%     if isempty(ind)  
%         kk = 0.5;
%         while isempty(ind)
%             nearby_lats = find(abs(slat-flat) < (kk+0.5)*target_res);
%             nearby_lons = find(abs(slon-flon) < (kk+0.5)*target_res);
%             ind = find(ismember(nearby_lats, nearby_lons));
%             kk = kk + 0.5;
%         end
%         ind = ind(1);
%         disp(['Filling in missing data for lon: ' num2str(slon(nearby_lons(ind))) ...
%             ' lat: ' num2str(slat(nearby_lats(ind))) ' with data from a nearby cell.'])
%         disp(['Final search radius = ' num2str((kk+0.5)*target_res) 'degrees'])
%         nmissing = nmissing + 1;     
%     end
        
    cellnumber = nearby_lats(ind);
    newlat = slat(cellnumber);
    newlon = slon(cellnumber);
    forcname_new = ['Forcings_' num2str(newlat, precision) '_' num2str(newlon, precision)];
    copyfile(fullfile(forcingdir, forcnames(k).name), fullfile(copyloc, forcname_new))
    
end

disp(['There were ' num2str(nmissing) 'grid cells of of: ' num2str(ncells) 'cells'])

return
