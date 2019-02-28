% Loads VIC4 or VIC5 classic vegetation parameter file into Matlab in a
% convenient, easy-to-work-with structure format
% Feb. 28, 2019 JRS
%
% Assumes that LAI is provided in the vegetation parameter file, but not
% albedo, etc.
%

% cd /Volumes/HD3/VICParametersGlobal/vic_params_global_0.5deg/
% vegfile = 'vegparam_test.txt';
% vegfile = 'global_veg_param_new';

function VEGPAR = load_veg_parameters(vegfile)

% load the vegetation parameter file
% maybe it would be useful to split it into two files?
% there might be a tonic script to convert vegetation parameter file to
% netcdf?

% initialize variables

% load names of vegetation types (get them from the vegetation library)
fID = fopen('vegnames.txt');
vegnames = textscan(fID, '%s');
vegnames = vegnames{1};
fclose(fID);

ncells = 61345;
VEGPAR(ncells) = struct(); % initializing a structure to a certain size

for i=1:ncells
    VEGPAR(i).cellID = [];
    VEGPAR(i).nveg = [];
    
    for j=1:length(vegnames)
        VEGPAR(i).(vegnames{j}).cv = [];
        VEGPAR(i).(vegnames{j}).rootdepth = zeros(0,2);
        VEGPAR(i).(vegnames{j}).rootfract = zeros(0,2);
        VEGPAR(i).(vegnames{j}).rootfract = zeros(0,2);
        VEGPAR(i).(vegnames{j}).LAI = zeros(0,12);
    end
    
end

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
      
        if current_cellID > 1
            VEGPAR(current_cellID).cellID = [VEGPAR(current_cellID).cellID, current_cellID]; 
            VEGPAR(current_cellID).nveg = [VEGPAR(current_cellID).nveg, current_nveg];  
        else
            VEGPAR(current_cellID).cellID = current_cellID;
            VEGPAR(current_cellID).nveg = current_nveg;  
        end
        
        prev_has_veg = 0;
        
    elseif nline(2) >= 1
%         disp(['There are ' num2str(nline(2)) ' vegetation classes in cell ' num2str(nline(1))])

        current_cellID = nline(1);
        current_nveg = nline(2);
      
        if current_cellID > 1
            VEGPAR(current_cellID).cellID = [VEGPAR(current_cellID).cellID, current_cellID]; 
            VEGPAR(current_cellID).nveg = [VEGPAR(current_cellID).nveg, current_nveg];  
        else
            VEGPAR(current_cellID).cellID = current_cellID;
            VEGPAR(current_cellID).nveg = current_nveg;  
        end
                
        tline = fgetl(fID); % read the next line, which has vegclass, cv, etc.
        if tline==-1; break; end
        
        nline = str2num(tline);
        currentlinenumber = currentlinenumber + 1;
        while length(nline)>2
%             disp('Reading data for each vegetation class in the grid cell')
            
            current_vegclass = nline(1); % vegclass numbers must correctly line up with the lines in vegnames
            current_cv = nline(2);
            current_rootdepth = [nline(3), nline(4)];
            current_rootfract = [nline(5), nline(6)];
            
            VEGPAR(current_cellID).(vegnames{current_vegclass}).cv = current_cv;
            VEGPAR(current_cellID).(vegnames{current_vegclass}).rootdepth = current_rootdepth;
            VEGPAR(current_cellID).(vegnames{current_vegclass}).rootfract = current_rootfract;
            
            tline = fgetl(fID); % read the next line, which has LAI
            if tline==-1; break; end
            
            nline = str2num(tline);
            currentlinenumber = currentlinenumber + 1;
            
            VEGPAR(current_cellID).(vegnames{current_vegclass}).LAI = ...
                [nline(1), nline(2), nline(3), nline(4), nline(5), nline(6), ...
                nline(7), nline(8), nline(9), nline(10), nline(11), nline(12)];
                        
            tline = fgetl(fID);
            if tline==-1; break; end
            
            nline = str2num(tline); 
            currentlinenumber = currentlinenumber + 1;
            prev_has_veg = 1;
            
        end        
        
    end
    
        if ~prev_has_veg % because then it has already advanced the line
            
            tline = fgetl(fID);
            if tline==-1; break; end
            
            nline = str2num(tline);
            currentlinenumber = currentlinenumber + 1;
        end
    
        disp(num2str(round(100*currentlinenumber/471929)))
        
end
fclose(fID);

% Save the processed vegetation parameter data
save('global_vegetation_parameters.mat', 'VEGPAR')
% It would be useful to be able to also detect which vegetation classes
% have values for a particular cell, but this is nonessential for now.

return 
