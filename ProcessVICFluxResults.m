function FLUXES = ProcessVICFluxResults(gridcells, results, nlayers, run_type, rec_interval)

% After VIC results are loaded into Matlab, as in LoadVICResults, this
% function processes the results, putting them into a structure with entries for each 
% grid cell. Adapted from Keith Cherkaur's 1998 S code.
%
% Note: Assumes time steps are the same for each flux output.
% Note: Only supports daily recording intervals
% Note: Only supports water balance run_type
%
% Todo: Somehow need to map the number of the gridcell to the lat/lon.
% 
% INPUTS
% gridcells = locations of the gridcells (lat/lon). Output from LoadVICResults.
% results (format should be an array with dimensions of n_timesteps x
% n_variables x n_gridcells)
% run_type = string (WATER_BALANCE, FULL_ENERGY, or FROZEN_SOIL)
%
% OUTPUTS
% Struct of VIC flux results.

ncells = size(results,3);

FLUXES.time = GetDateTime(results);

for k=1:ncells

    if strcmp(run_type, 'WATER_BALANCE')
        
        if strcmp(rec_interval, 'hourly')
            FLUXES = NaN;
                    
        elseif strcmp(rec_interval, 'daily')
            
            prec = results(:,4,k);
            evap = results(:,5,k);
            runoff = results(:,6,k);
            baseflow = results(:,7,k);
            wdew = results(:,8,k);
            moist = results(:,9:nlayers + 8,k);
            net_short = results(:,nlayers + 9,k);
            r_net = results(:,nlayers + 10,k);
            evap_canop = results(:,nlayers + 11,k);
            evap_veg = results(:,nlayers + 12,k);
            evap_bare = results(:,nlayers + 13,k);
            sub_canop = results(:,nlayers + 14,k);
            sub_snow = results(:,nlayers + 15,k);
            aero_resist = results(:,nlayers + 16,k);
            surf_temp = results(:,nlayers + 17,k);
            albedo = results(:,nlayers + 18,k);
            rel_humid = results(:,nlayers + 19,k);
            in_long = results(:,nlayers + 20,k);
            air_temp = results(:,nlayers + 21,k);
            wind = results(:,nlayers + 22,k);
            
            T = table(prec, evap, runoff, baseflow, wdew, moist, net_short, ...
                r_net, evap_canop, evap_veg, evap_bare, sub_canop, sub_snow, ...
                aero_resist, surf_temp, albedo, rel_humid, in_long, air_temp, wind);
            % tmp = gridcells(k).name;
            cmd = ['FLUXES.gridcell_' gridcells{k} ' = T;'];
            eval(cmd);
            
        elseif strcmp(rec_interval, 'monthly')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'yearly')
            FLUXES = NaN;
        end
        
    elseif strcmp(run_type, 'FULL_ENERGY')
        
        if strcmp(rec_interval, 'hourly')
            FLUXES = NaN; 
        elseif strcmp(rec_interval, 'daily')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'monthly')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'yearly')
            FLUXES = NaN;
        end
        
    elseif strcmp(run_type, 'FROZEN_SOIL')
        if strcmp(rec_interval, 'hourly')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'daily')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'monthly')
            FLUXES = NaN;
        elseif strcmp(rec_interval, 'yearly')
            FLUXES = NaN;
        end
        
    end

end
    
if ~isstruct(FLUXES)
    warning('ProcessVICFluxResults only supports daily flux outputs at this time')
    warning('ProcessVICFluxResults only supports WATER_BALANCE run_type at this time')
else
    FLUXES.units.prec =  'mm';
    FLUXES.units.evap =  'mm';
    FLUXES.units.runoff =  'mm';
    FLUXES.units.baseflow =  'mm';
    FLUXES.units.wdew =  'mm';
    FLUXES.units.moist =  'mm';
    FLUXES.units.net_short =  'W/m^2';
    FLUXES.units.r_net =  'W/m^2';
    FLUXES.units.evap_canop =  'mm';
    FLUXES.units.evap_veg =  'mm';
    FLUXES.units.evap_bare =  'mm';
    FLUXES.units.sub_canop =  'mm';
    FLUXES.units.sub_snow =  'mm';
    FLUXES.units.aero_resist =  's/m';
    FLUXES.units.surf_temp =  'deg. C';
    FLUXES.units.albedo =  '-';
    FLUXES.units.rel_humid =  '-';
    FLUXES.units.in_long =  'W/m^2';
    FLUXES.units.air_temp = 'deg. C';
    FLUXES.units.wind = 'm/s';
end

end