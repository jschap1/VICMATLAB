function SNOW = ProcessVICSnowResults(gridcells, results, run_type, rec_interval)

% After VIC results are loaded into Matlab, as in LoadVICResults, this
% function processes the results, putting them into a structure with entries for each 
% grid cell. Similar to ProcessVICFluxResults.
%
% INPUTS
% gridcells = locations of the gridcells (lat/lon). Output from LoadVICResults.
% results (format should be an array with dimensions of n_timesteps x
% n_variables x n_gridcells)
% run_type = string (WATER_BALANCE, FULL_ENERGY, or FROZEN_SOIL)
%
% OUTPUTS
% Struct of VIC snow results.

ncells = size(results,3);

for k=1:ncells
%% water balance
    if strcmp(run_type, 'WATER_BALANCE') % this is a more compact (better) format than ProcessVICFluxResults
         
       if strcmp(rec_interval, 'hourly')   
            Y = results(:,1,1);
            M = results(:,2,1);
            D = results(:,3,1);
            H = results(:,4,1);
            SNOW.time = datetime(Y,M,D,H,0,0);  
            t_prec = 4; % index for date precision
        elseif strcmp(rec_interval, 'daily') 
            Y = results(:,1,1);
            M = results(:,2,1);
            D = results(:,3,1);
            SNOW.time = datetime(Y,M,D);  
            t_prec = 3;
        elseif strcmp(rec_interval, 'monthly')
            Y = results(:,1,1);
            M = results(:,2,1);
            SNOW.time = datetime(Y,M,0);  
            t_prec = 2;
        elseif strcmp(rec_interval, 'yearly')
            Y = results(:,1,1);
            SNOW.time = datetime(Y,0,0);  
            t_prec = 1;
        end
        
        swe = results(:,t_prec + 1,k);
        snow_depth = results(:,t_prec + 2,k);
        snow_canopy = results(:,t_prec + 3,k);
        snow_cover = results(:,t_prec + 4,k);

        T = table(swe, snow_depth, snow_canopy, snow_cover);          
        SNOW.ts.(gridcells{k}) = T;   

%% full energy
    elseif strcmp(run_type, 'FULL_ENERGY')
        
        if strcmp(rec_interval, 'hourly')
            SNOW = NaN; 
        elseif strcmp(rec_interval, 'daily')
            SNOW = NaN;
        elseif strcmp(rec_interval, 'monthly')
            SNOW = NaN;
        elseif strcmp(rec_interval, 'yearly')
            SNOW = NaN;
        end
     
%% frozen soil      
    elseif strcmp(run_type, 'FROZEN_SOIL')
        if strcmp(rec_interval, 'hourly')
            SNOW = NaN;
        elseif strcmp(rec_interval, 'daily')
            SNOW = NaN;
        elseif strcmp(rec_interval, 'monthly')
            SNOW = NaN;
        elseif strcmp(rec_interval, 'yearly')
            SNOW = NaN;
        end
        
    end

end
    
%%
if ~isstruct(SNOW)
    warning('ProcessVICSnowResults only supports WATER_BALANCE run_type at this time')
else
    SNOW.units.swe =  'mm';
    SNOW.units.snow_depth =  'cm';
    SNOW.units.snow_canopy =  'mm';
    SNOW.units.snow_cover =  '-';
end

return