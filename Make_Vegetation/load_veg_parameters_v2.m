% Loads VIC4 or VIC5 classic vegetation parameter file
% May 2, 2019 JRS
%
% Loads vegetation parameter file into Matlab as a relational database
% Assumes that LAI is provided in the vegetation parameter file, but not albedo, etc.
% 
% INPUTS
% vegnamefile, names of vegetation types (get them from the vegetation library)
% soilfile, soil parameter file
% vegfile, vegetation parameter file
%
% OUTPUTS
% nvegtable, table with grid cell IDs and nveg
% vegparamtable, table with cell IDs, cover types, and vegetation
% parameters
% latlontable, table with cell IDs and lat/lons
%
% SAMPLE INPUTS
% soilfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_1/soils_3L_MERIT.txt';
% vegfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_1/global_vegetation_1_16.txt';

function [nvegtable, vegparamtable, latlontable] = load_veg_parameters_v2(vegfile, soilfile)

% get lat/lon and ncells
soils = load(soilfile);
soils = soils(:,1:5);
ncells = size(soils, 1);
latlontable = [soils(:,2) soils(:,3:4)];
disp('Read soil parameter file')

%% Read through the whole vegetation parameter file, line by line

nvegtable = zeros(ncells, 2); % number of vegetation classes per grid cell

fID = fopen(vegfile);
tline = fgetl(fID);
nline = str2num(tline);
currentlinenumber = 1; % keep track of line number
prev_has_veg = 0; % for indexing purposes
while ischar(tline)
    
    if nline(2) == 0
%         disp(['No vegetation classes in cell ' num2str(nline(1))])

        current_cellID = nline(1);
        current_nveg = nline(2);
        nvegtable(current_cellID, 1) = current_cellID;
        nvegtable(current_cellID, 2) = current_nveg;  
        prev_has_veg = 0;
        
    elseif nline(2) >= 1
%         disp(['There are ' num2str(nline(2)) ' vegetation classes in cell ' num2str(nline(1))])

        current_cellID = nline(1);
        current_nveg = nline(2);
      
        nvegtable(current_cellID, 1) = current_cellID; 
        nvegtable(current_cellID, 2) = current_nveg;  
    
        tline = fgetl(fID); % read the next line, which has vegclass, cv, etc.
        if ~ischar(tline); break; end
                
        nline = str2num(tline);
        currentlinenumber = currentlinenumber + 1;
        while length(nline)>2
%             disp('Reading data for each vegetation class in the grid cell')
                        
            tline = fgetl(fID); % read the next line, which has LAI
            if ~ischar(tline); break; end
            
            nline = str2num(tline);
            currentlinenumber = currentlinenumber + 1;
                                    
            tline = fgetl(fID);
            if ~ischar(tline); break; end
            
            nline = str2num(tline); 
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
        disp(['progress: ' num2str(round(100*current_cellID/ncells)) '%'])
    end
         
end
fclose(fID);

disp('Read through vegetation parameter file once')
disp('One more time is necessary to extract all data')

%% Initialize variables

vegparamtable = zeros(sum(nvegtable(:,2)), 19);

% Make header for vegparamtable
LAI_string = cell(12, 1);
for k=1:12
    LAI_string{k} = ['LAI' num2str(k)];
end
other_header_names = {'cover_type','cover_fraction','rootfract1','rootfract2','rootdepth1','rootdepth2'};
header_names = cell(19, 1);
for k=1:19
    if k > 7
        header_names{k} = LAI_string{k-7};
    elseif k==1
        header_names{k} = 'cell_ID';
    else
        header_names{k} = other_header_names{k-1};
    end
end

%% Read through the vegetation parameter file again

fID = fopen(vegfile);
tline = fgetl(fID);
vpt_index = 1;
currentlinenumber = 1;

while ischar(tline)
    
    nline = str2num(tline); % this takes the most computation time
    current_cellID = nline(1);
    
    if nline(2) >= 1
        % has vegetation
        vegparamtable(vpt_index, 1) = current_cellID;
        tline = fgetl(fID);
        nline = str2num(tline);
        currentlinenumber = currentlinenumber + 1;
        while length(nline)>2
%             disp('Reading data for each vegetation class in the grid cell')
            
            current_vegclass = nline(1); % vegclass numbers must correctly line up with the lines in vegnames
            current_cv = nline(2);
            current_rootdepth = [nline(3), nline(4)];
            current_rootfract = [nline(5), nline(6)];
            
            vegparamtable(vpt_index, 2) = current_vegclass; % cover type
            vegparamtable(vpt_index, 3) = current_cv; % cover fraction
            vegparamtable(vpt_index, 4:5) = current_rootfract;
            vegparamtable(vpt_index, 6:7) = current_rootdepth;
             
            tline = fgetl(fID); % read the next line, which has LAI           
            nline = str2num(tline);
            currentlinenumber = currentlinenumber + 1;
            
            vegparamtable(vpt_index, 8:19) = [nline(1), nline(2), ...
                nline(3), nline(4), nline(5), nline(6), ...
                nline(7), nline(8), nline(9), nline(10), nline(11), nline(12)];
           
            tline = fgetl(fID);
            nline = str2num(tline); 
            currentlinenumber = currentlinenumber + 1;            
            vpt_index = vpt_index + 1; % update vpt row index
        end          
    else
        % does not have vegetation
        tline = fgetl(fID);
        currentlinenumber = currentlinenumber + 1;
    end
    
    if mod(current_cellID, 1e4)==0
        disp(['progress: ' num2str(round(100*current_cellID/ncells)) '%'])
    end    
    
end

%% Save results 

save('./global_vegetation_parameters.mat', 'nvegtable', 'vegparamtable', 'latlontable')
disp(['Saved ' 'global_vegetation_parameters.mat'])

% It would be useful to be able to also detect which vegetation classes
% have values for a particular cell, but this is nonessential for now.

return 
