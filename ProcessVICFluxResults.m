function FLUXES = ProcessVICFluxResults(gridcells, results, nlayers, run_type, rec_interval)

% After VIC results are loaded into Matlab, as in LoadVICResults, this
% function processes the results, putting them into a structure with entries for each 
% grid cell. Adapted from Keith Cherkaur's 1998 S code.
%
% Note: Assumes time steps are the same for each flux output.
% Note: Only supports daily recording intervals
% Note: Only supports water balance run_type
%
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

for k=1:ncells
%% water balance
    if strcmp(run_type, 'WATER_BALANCE')
        
        if strcmp(rec_interval, 'hourly')
            
            Y = results(:,1,1);
            M = results(:,2,1);
            D = results(:,3,1);
            H = results(:,4,1);
            FLUXES.time = datetime(Y,M,D,H,0,0);         
            
            prec = results(:,5,k);
            evap = results(:,6,k);
            runoff = results(:,7,k);
            baseflow = results(:,8,k);
            wdew = results(:,9,k);
            moist = results(:,10:nlayers + 9,k);
            net_short = results(:,nlayers + 10,k);
            r_net = results(:,nlayers + 11,k);
            evap_canop = results(:,nlayers + 12,k);
            evap_veg = results(:,nlayers + 13,k);
            evap_bare = results(:,nlayers + 14,k);
            sub_canop = results(:,nlayers + 15,k);
            sub_snow = results(:,nlayers + 16,k);
            aero_resist = results(:,nlayers + 17,k);
            surf_temp = results(:,nlayers + 18,k);
            albedo = results(:,nlayers + 19,k);
            rel_humid = results(:,nlayers + 20,k);
            in_long = results(:,nlayers + 21,k);
            air_temp = results(:,nlayers + 22,k);
            wind = results(:,nlayers + 23,k);
            
            T = table(prec, evap, runoff, baseflow, wdew, moist, net_short, ...
                r_net, evap_canop, evap_veg, evap_bare, sub_canop, sub_snow, ...
                aero_resist, surf_temp, albedo, rel_humid, in_long, air_temp, wind);          
            FLUXES.ts.(gridcells{k}) = T;
                    
        elseif strcmp(rec_interval, 'daily')
            
            Y = results(:,1,1);
            M = results(:,2,1);
            D = results(:,3,1);
            FLUXES.time = datetime(Y,M,D); 
            
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
            FLUXES.ts.(gridcells{k}) = T;
            
        elseif strcmp(rec_interval, 'monthly')

            Y = results(:,1,1);
            M = results(:,2,1);
            FLUXES.time = datetime(Y,M,0);         
            
            prec = results(:,3,k);
            evap = results(:,4,k);
            runoff = results(:,5,k);
            baseflow = results(:,6,k);
            wdew = results(:,7,k);
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
            FLUXES.ts.(gridcells{k}) = T;            
            
        elseif strcmp(rec_interval, 'yearly')
            
            Y = results(:,1,1);
            FLUXES.time = datetime(Y,0,0);         
            
            prec = results(:,2,k);
            evap = results(:,3,k);
            runoff = results(:,4,k);
            baseflow = results(:,5,k);
            wdew = results(:,6,k);
            moist = results(:,8:nlayers + 7,k);
            net_short = results(:,nlayers + 8,k);
            r_net = results(:,nlayers + 9,k);
            evap_canop = results(:,nlayers + 10,k);
            evap_veg = results(:,nlayers + 11,k);
            evap_bare = results(:,nlayers + 12,k);
            sub_canop = results(:,nlayers + 13,k);
            sub_snow = results(:,nlayers + 14,k);
            aero_resist = results(:,nlayers + 15,k);
            surf_temp = results(:,nlayers + 16,k);
            albedo = results(:,nlayers + 17,k);
            rel_humid = results(:,nlayers + 18,k);
            in_long = results(:,nlayers + 19,k);
            air_temp = results(:,nlayers + 20,k);
            wind = results(:,nlayers + 21,k);
            
            T = table(prec, evap, runoff, baseflow, wdew, moist, net_short, ...
                r_net, evap_canop, evap_veg, evap_bare, sub_canop, sub_snow, ...
                aero_resist, surf_temp, albedo, rel_humid, in_long, air_temp, wind);          
            FLUXES.ts.(gridcells{k}) = T;               
            
        end

%% full energy
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
     
%% frozen soil      
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
    
%%
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