
%%
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

% cd /Volumes/HD3/SWOTDA/Data/IRB/VIC
% forcnames = dir('./Forc_halfyear/Forcings_*');
% soils = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/soils.SB');
% 
% slat = soils(:,3);
% slon = soils(:,4);

% find closest
target_res = 1/16;
precision = '%3.5f';
copyloc = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/Forc_halfyear_input';

for k=1:length(forcnames)
    
    tmp = strsplit(forcnames(k).name, '_');
    flat = str2double(tmp{2});
    tmp = strsplit(tmp{3}, '.txt');
    flon = str2double(tmp{1});

    nearby_lats = find(abs(slat-flat) < 0.5*target_res);
    nearby_lons = find(abs(slon-flon) < 0.5*target_res);
   
    ind = find(ismember(nearby_lats, nearby_lons));
    
    if isempty(ind)
        continue
    end
    
    cellnumber = nearby_lats(ind);
    newlat = slat(cellnumber);
    newlon = slon(cellnumber);
    
    forcname_new = ['Forcings_' num2str(newlat, precision) '_' num2str(newlon, precision)];
    copyfile(fullfile('./Forc_halfyear/', forcnames(k).name), fullfile(copyloc, forcname_new))
     
    disp(forcname_new)

    % report progress
%     if mod(k, 1e3) == 0
%         disp(k)
%     end
    
end
