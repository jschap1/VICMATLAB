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
    FLUXES.time = datetime(Y,M,D,H,0,0);   
    case 'daily'
    t_prec = 3;
    Y = results(:,1,1);
    M = results(:,2,1);
    D = results(:,3,1);
    FLUXES.time = datetime(Y,M,D);  
    case 'monthly'
    t_prec = 2; 
    Y = results(:,1,1);
    M = results(:,2,1);
    FLUXES.time = datetime(Y,M,0);  
    case 'yearly'
    t_prec = 1;
    Y = results(:,1,1);
    FLUXES.time = datetime(Y,0,0);  
    otherwise
        error('Specify recording interval')
end
%%%%%%

for k=1:ncells
%% water balance
    if strcmp(run_type, 'WATER_BALANCE')
                       
        prec = results(:,t_prec + 1,k);
        evap = results(:,t_prec + 2,k);
        runoff = results(:,t_prec + 3,k);
        baseflow = results(:,t_prec + 4,k);
        wdew = results(:,t_prec + 5,k);
        moist = results(:,t_prec + 6:nlayers + t_prec + 5,k);
        net_short = results(:,nlayers + t_prec + 6,k);
        r_net = results(:,nlayers + t_prec + 7,k);
        evap_canop = results(:,nlayers + t_prec + 8,k);
        evap_veg = results(:,nlayers + t_prec + 9,k);
        evap_bare = results(:,nlayers + t_prec + 10,k);
        sub_canop = results(:,nlayers + t_prec + 11,k);
        sub_snow = results(:,nlayers + t_prec + 12,k);
        aero_resist = results(:,nlayers + t_prec + 13,k);
        surf_temp = results(:,nlayers + t_prec + 14,k);
        albedo = results(:,nlayers + t_prec + 15,k);
        rel_humid = results(:,nlayers + t_prec + 16,k);
        in_long = results(:,nlayers + t_prec + 17,k);
        air_temp = results(:,nlayers + t_prec + 18,k);
        wind = results(:,nlayers + t_prec + 19,k);

        T = table(prec, evap, runoff, baseflow, wdew, moist, net_short, ...
            r_net, evap_canop, evap_veg, evap_bare, sub_canop, sub_snow, ...
            aero_resist, surf_temp, albedo, rel_humid, in_long, air_temp, wind);          
        FLUXES.ts.(gridcells{k}) = T;

%% full energy
    elseif strcmp(run_type, 'FULL_ENERGY')
        
        prec = results(:,t_prec + 1,k);
        evap = results(:,t_prec + 2,k);
        runoff = results(:,t_prec + 3,k);
        baseflow = results(:,t_prec + 4,k);
        wdew = results(:,t_prec + 5,k);
        moist = results(:,t_prec + 6:nlayers + t_prec + 5,k);
        rad_temp = results(:,nlayers + t_prec + 6,k); %
        net_short = results(:,nlayers + t_prec + 7,k);
        r_net = results(:,nlayers + t_prec + 8,k);
        latent = results(:,nlayers + t_prec + 9,k); %
        evap_canop = results(:,nlayers + t_prec + 10,k);
        evap_veg = results(:,nlayers + t_prec + 11,k);
        evap_bare = results(:,nlayers + t_prec + 12,k);
        sub_canop = results(:,nlayers + t_prec + 13,k);
        sub_snow = results(:,nlayers + t_prec + 14,k);
        sensible = results(:,nlayers + t_prec + 15,k); %
        grnd_flux = results(:,nlayers + t_prec + 16,k); %
        deltah = results(:,nlayers + t_prec + 17,k); %
        fusion = results(:,nlayers + t_prec + 18,k); %
        aero_resist = results(:,nlayers + t_prec + 19,k);
        surf_temp = results(:,nlayers + t_prec + 20,k);
        albedo = results(:,nlayers + t_prec + 21,k);
        rel_humid = results(:,nlayers + t_prec + 22,k);
        in_long = results(:,nlayers + t_prec + 23,k);
        air_temp = results(:,nlayers + t_prec + 24,k);
        wind = results(:,nlayers + t_prec + 25,k);

        T = table(prec, evap, runoff, baseflow, wdew, moist, rad_temp, ...
            net_short, r_net, latent, evap_canop, evap_veg, evap_bare, ...
            sub_canop, sub_snow, sensible, grnd_flux, deltah, fusion, ...
            aero_resist, surf_temp, albedo, rel_humid, in_long, ...
            air_temp, wind);          
        FLUXES.ts.(gridcells{k}) = T;        
     
%% frozen soil      
    elseif strcmp(run_type, 'FROZEN_SOIL')
        FLUXES = NaN;
        
    end

end
    
%%
if ~isstruct(FLUXES)
    warning('ProcessVICFluxResults only supports daily flux outputs at this time')
    warning('ProcessVICFluxResults only supports WATER_BALANCE run_type at this time')
else % define units. These are the default VIC units, not ALMA units
    FLUXES.units.prec =  'mm';
    FLUXES.units.evap =  'mm';
    FLUXES.units.runoff =  'mm';
    FLUXES.units.baseflow =  'mm';
    FLUXES.units.wdew =  'mm';
    FLUXES.units.moist =  'mm';
    FLUXES.units.rad_temp = 'K';
    FLUXES.units.net_short =  'W/m^2';
    FLUXES.units.r_net =  'W/m^2';
    FLUXES.units.latent = 'W/m^2'; 
    FLUXES.units.evap_canop =  'mm';
    FLUXES.units.evap_veg =  'mm';
    FLUXES.units.evap_bare =  'mm';
    FLUXES.units.sub_canop =  'mm';
    FLUXES.units.sub_snow =  'mm';
    FLUXES.units.sensible = 'W/m^2'; 
    FLUXES.units.grnd_flux = 'W/m^2'; 
    FLUXES.units.deltah = 'W/m^2'; 
    FLUXES.units.fusion = 'W/m^2'; 
    FLUXES.units.aero_resist =  's/m';
    FLUXES.units.surf_temp =  'deg. C';
    FLUXES.units.albedo =  '-';
    FLUXES.units.rel_humid =  '-';
    FLUXES.units.in_long =  'W/m^2';
    FLUXES.units.air_temp = 'deg. C';
    FLUXES.units.wind = 'm/s';
end

end