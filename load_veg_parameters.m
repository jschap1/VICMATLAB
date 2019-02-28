% Loads VIC4 or VIC5 classic vegetation parameter file into Matlab in a
% convenient, easy-to-work-with structure format

% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % no organic, no frozen soil, no July_Tavg, two soil layers
%     'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
%     'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
%     'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
%     'soil_dens1','soil_dens2', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2'};

% writecsv_cell('varnames.txt', varnames);

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

% cellID = [];
% nveg = [];
% vegclass = [];
% cv = [];
% rootdepth1 = [];
% rootdepth2 = [];
% rootfract1 = [];
% rootfract2 = [];
% LAI1 = [];
% LAI2 = [];
% LAI3 = [];
% LAI4 = [];
% LAI5 = [];
% LAI6 = [];
% LAI7 = [];
% LAI8 = [];
% LAI9 = [];
% LAI10 = [];
% LAI11 = [];
% LAI12 = [];



fID = fopen(vegfile);
tline = fgetl(fID);
nline = str2num(tline);
currentlinenumber = 1; % keep track of line number
prev_has_veg = 0; % for indexing purposes
while ischar(tline)
    
    if nline(2) == 0
%         disp(['No vegetation classes in cell ' num2str(nline(1))])

%         cellID = [cellID, nline(1)];
%         nveg = [nveg, nline(2)]; 

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

%         cellID = [cellID, nline(1)];
%         nveg = [nveg, nline(2)];

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
            
%             vegclass = [vegclass nline(1)];
%             cv = [cv nline(2)];
%             rootdepth1 = [rootdepth1 nline(3)];
%             rootdepth2 = [rootdepth2 nline(4)];
%             rootfract1 = [rootfract1 nline(5)];
%             rootfract2 = [rootfract2 nline(6)];

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
            
%             LAI1 = [LAI1 nline(1)];
%             LAI2 = [LAI2 nline(2)];
%             LAI3 = [LAI3 nline(3)];
%             LAI4 = [LAI4 nline(4)];
%             LAI5 = [LAI5 nline(5)];
%             LAI6 = [LAI6 nline(6)];
%             LAI7 = [LAI7 nline(7)];
%             LAI8 = [LAI8 nline(8)];
%             LAI9 = [LAI9 nline(9)];
%             LAI10 = [LAI10 nline(10)];
%             LAI11 = [LAI11 nline(11)];
%             LAI12 = [LAI12 nline(12)];
            
            % check if next line is a new tile or not
            
%             C = fseek(fID, 1, 1)
            
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
