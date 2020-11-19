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

%%%%%%
% Assumes timesteps are the same for each grid cell
% Otherwise, move this inside the k loop.
switch rec_interval 
    case 'hourly'
    t_prec = 4; % index for date precision
    Y = results(:,1,1);
    M = results(:,2,1);
    D = results(:,3,1);
    H = results(:,4,1);
    SNOW.time = datetime(Y,M,D,H,0,0);   
    case 'daily'
    t_prec = 3;
    Y = results(:,1,1);
    M = results(:,2,1);
    D = results(:,3,1);
    SNOW.time = datetime(Y,M,D);  
    case 'monthly'
    t_prec = 2; 
    Y = results(:,1,1);
    M = results(:,2,1);
    SNOW.time = datetime(Y,M,0);  
    case 'yearly'
    t_prec = 1;
    Y = results(:,1,1);
    SNOW.time = datetime(Y,0,0);  
    otherwise
        error('Specify recording interval')
end
%%%%%%

for k=1:ncells
%% water balance
    if strcmp(run_type, 'WATER_BALANCE')
                
        swe = results(:,t_prec + 1,k);
        snow_depth = results(:,t_prec + 2,k);
        snow_canopy = results(:,t_prec + 3,k);
        snow_cover = results(:,t_prec + 4,k);

        T = table(swe, snow_depth, snow_canopy, snow_cover);          
        SNOW.ts.(gridcells{k}) = T;   

%% full energy
    elseif strcmp(run_type, 'FULL_ENERGY')
        
        swe = results(:,t_prec + 1,k);
        snow_depth = results(:,t_prec + 2,k);
        snow_canopy = results(:,t_prec + 3,k);
        snow_cover = results(:,t_prec + 4,k);
        
        advection = results(:,t_prec + 5,k);
        deltacc = results(:,t_prec + 6,k);
        snow_flux = results(:,t_prec + 7,k);
        rfrz_energy = results(:,t_prec + 8,k);
        melt_energy = results(:,t_prec + 9,k);
        adv_sens = results(:,t_prec + 10,k);
        latent_sub = results(:,t_prec + 11,k);
        snow_surf_temp = results(:,t_prec + 12,k);
        snow_pack_temp = results(:,t_prec + 13,k);
        snow_melt = results(:,t_prec + 14,k);
        
        T = table(swe, snow_depth, snow_canopy, snow_cover, advection, ...
            deltacc, snow_flux, rfrz_energy, melt_energy, adv_sens, ...
            latent_sub, snow_surf_temp, snow_pack_temp, snow_melt);          
        SNOW.ts.(gridcells{k}) = T;
     
%% frozen soil      
    elseif strcmp(run_type, 'FROZEN_SOIL')
        SNOW = NaN;
        
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
    SNOW.units.advection = 'W/m^2';
    SNOW.units.deltacc = 'W/m^2';
    SNOW.units.snow_flux = 'W/m^2';
    SNOW.units.rfrz_energy = 'W/m^2';
    SNOW.units.melt_energy = 'W/m^2';
    SNOW.units.adv_sens = 'W/m^2';         
    SNOW.units.latent_sub = 'W/m^2';
    SNOW.units.snow_surf_temp = 'deg C';
    SNOW.units.snow_pack_temp = 'deg C'; 
    SNOW.units.snow_melt = 'mm';
end

return