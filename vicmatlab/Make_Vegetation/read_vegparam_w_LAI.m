% Read vegetation parameter file
%
% Loads the VIC4 or VIC5 classic vegetation parameter file
% 2/18/2020 JRS
% Modified from read_vegparam to accommodate monthly LAI 
%
% Loads vegetation parameter file into Matlab as a relational database
%
% INPUTS
% soilfile, soil parameter file
% vegfile, vegetation parameter file
% key, names of vegetation types (get them from the vegetation library)
%
% OUTPUTS
% nvegtable, table with grid cell IDs and nveg
% vegparamtable, table with cell IDs, cover types, and vegetation
% parameters
% nlayers = number of soil layers
% latlontable, table with cell IDs and lat/lons
% ncells = number of cells in the vegetation parameter file (if known)
%
% SAMPLE INPUTS
% soilfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_3/soils_3L_MERIT.txt';
% vegfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_3/global_vegetation_1_16_sub.txt';
%
% The user has the option to input soil data directly, instead of the name of the soil
% parameter file, to avoid having to load the soil parameter file if it is
% already loaded

function [nvegtable, vegparamtable, latlontable, LC] = read_vegparam_w_LAI(vegfile, soils, ncells, savename)

if ischar(soils)
    soils = load(soils);
    disp('Read soil parameter file')
end

% get lat/lon and ncells
soils = soils(:,1:5);
latlontable = [soils(:,2) soils(:,3:4)];
% ncells = size(soils, 1);

%% Read through the whole vegetation parameter file, line by line
% The first read-through is to obtain the number of vegetation classes for
% each grid cell (and the number of cells in the vegetation parameter file)

if isempty(ncells)
    nvegtable = zeros(0, 2); % number of vegetation classes per grid cell
    flag = 0;
else
    nvegtable = zeros(ncells, 2);
    flag = 1; % flag for whether or not ncells is known a priori
end

k = 0; % index for grid cell

fID = fopen(vegfile);
tline = fgetl(fID);
nline = str2num(tline);
currentlinenumber = 1; % keep track of line number
prev_has_veg = 0; % for indexing purposes
% nvegtable_ind = 1;

while ischar(tline)
    
    if nline(2) == 0
%         disp(['No vegetation classes in cell ' num2str(nline(1))])

        current_cellID = nline(1);
        current_nveg = nline(2);
        k = k + 1;
        
        nvegtable(k, 1) = current_cellID;
        nvegtable(k, 2) = current_nveg;  
        prev_has_veg = 0;
        
    elseif nline(2) >= 1
%         disp(['There are ' num2str(nline(2)) ' vegetation classes for cellID = ' num2str(nline(1))])

        current_cellID = nline(1);
        current_nveg = nline(2);
        k = k + 1;
        
        if flag
            nvegtable(k,:) = [current_cellID, current_nveg];
        else
            nvegtable = vertcat(nvegtable, [current_cellID, current_nveg]);
        end
        
        tline = fgetl(fID); % read the next line, which has vegclass, cv, etc.
        if ~ischar(tline); break; end
                
        nline = str2num(tline);
        currentlinenumber = currentlinenumber + 1;
        
        while length(nline)>2
%             disp('Reading data for each vegetation class in the grid cell')
       
            tline = fgetl(fID); % read the next line
            if ~ischar(tline); break; end            
            nline = str2num(tline); 
            if length(nline)==12
%                 disp('This line has LAI')
                  1;
            else
%                 disp('This line has vegetation parameters')
                  2;
            end
            currentlinenumber = currentlinenumber + 1;
            prev_has_veg = 1;
            
        end        
        
    end
    
    if ~prev_has_veg % because then it has already advanced the line

        tline = fgetl(fID);
        if ~ischar(tline); break; end
        nline = str2num(tline);
        currentlinenumber = currentlinenumber + 1;
        
    end

    if mod(current_cellID, 1e4)==0
        disp(['current_cellID: ' num2str(current_cellID)])
    end
             
end
fclose(fID);

if ~flag
    ncells = k-1; % number of grid cells in the vegetation parameter file
end

disp('Read through vegetation parameter file once')
disp('One more time is necessary to extract all data')

%% Initialize variables

% Initialize vegparamtable [cell ID, cover fraction, rf, rd] for each
% land cover class
% header_names = {'cover_fraction','rootfract1','rootfract2','rootfract3', 'rootdepth1','rootdepth2', 'rootdepth3'};
header_names = {'cover_fraction','rootfract1','rootfract2','rootfract3', 'rootdepth1','rootdepth2', 'rootdepth3', ...
    'LAI1','LAI2','LAI3','LAI4','LAI5','LAI6','LAI7','LAI8','LAI9','LAI10','LAI11','LAI12'};
nvars = length(header_names);

% classification_scheme = 'IGBP';
classification_scheme = 'UMD';
if strcmp(classification_scheme, 'UMD')
    
%     LC.class_names = {'Water_Bodies','Evergreen_Needleleaf', ...
%         'Evergreen_Broadleaf','Deciduous_Needleleaf','Deciduous_Broadleaf',...
%         'Mixed_Forests','Closed_Shrublands','Open_Shrublands','Savannas','Woody_Savannas',...
%         'Grasslands','Permanent_Wetlands','Croplands','Urban_and_Built_Up_Lands', ... 
%         'Cropland_Natural_Vegetation_Mosaics','Snow_and_Ice','Barren'};
    LC.class_names = {'Evergreen_Needleleaf','Evergreen_Broadleaf', ...
        'Deciduous_Needleleaf','Deciduous_Broadleaf',...
        'Mixed_Cover','Woodlands','Wooded_Grasslands','Closed_Shrublands','Open_Shrublands',...
        'Grasslands','Cropland'};
    LC.class_names = LC.class_names';
    LC.nclasses = length(LC.class_names);

%     LC.class_number = [17, 7, 13, 15, 5, 4, 3, 2, 11, 6, 8, 1, 12, 9, 16, 14, 10];
%     LC.alphabetical_order = [16 8 7 6 5 10 2 11 17 13 9 12 4 15 3 14 1];
%     LC.class_names = LC.class_names(LC.alphabetical_order)';
end

for i=1:LC.nclasses
    vegparamtable.(LC.class_names{i}) = zeros(ncells, nvars + 1);
end

%% Read through the vegetation parameter file again

% Re-writing the code to allow it to read in the Livneh vegetation
% parameters with LAI. Goal is to compare ET between VICGlobal and L2013
% simulations.
%
% vvvvvv--2/18/2020--vvvvvvvv


% This time through, populate the vegparamtable

fID = fopen(vegfile, 'r');
currentlinenumber = 0;
k = 0; % index for grid cell

while k<ncells
    
    tline = fgetl(fID);
    k = k + 1;
    currentlinenumber = currentlinenumber + 1;
    
    for i=1:LC.nclasses
        vegparamtable.(LC.class_names{i})(k, 1) = nvegtable(k,1); % cell ID
    end
    
    current_nveg = nvegtable(k, 2);
    for ii=1:current_nveg
        
        % Vegetation parameters
        tline = fgetl(fID);
        currentlinenumber = currentlinenumber + 1;
        nline = str2num(tline);
        current_vegtype = nline(1); % vegetation type number (IGBP code) 
        vegparamtable.(LC.class_names{current_vegtype})(k, 2:8) = nline(2:end);
        
        % LAI
        tline = fgetl(fID);
        currentlinenumber = currentlinenumber + 1;
        nline = str2num(tline);
        vegparamtable.(LC.class_names{current_vegtype})(k, 9:end) = nline(1:end); 
        
    end
    
    if mod(k, 1e4)==0
        disp(['progress: ' num2str(round(100*k/ncells)) '%'])
    end     
    
end

%% Save results 

save(savename, 'nvegtable', 'vegparamtable', 'latlontable', 'LC', '-v7.3')
disp(['Saved ' savename])

return 
