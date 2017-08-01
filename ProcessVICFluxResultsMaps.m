function FLUXES = ProcessVICFluxResultsMaps(FLUXES, t)

% Adds maps to the FLUXES structure generated in ProcessVICFluxResults.
%
% INPUTS
% Struct of VIC flux results, with time series only
% t = time index when you would like to make the map
%
% OUTPUTS
% Struct of VIC flux results, with maps and time series

ncells = length(fieldnames(FLUXES.ts));

    % Indexing convention
    % i, x
    % j, y
    % t, time
    % p, varnames
    % k, ncells

varnames = FLUXES.ts.fluxes_48_1875_120_6875.Properties.VariableNames;
cellnames = fieldnames(FLUXES.ts);

for p = 1:length(varnames)
    
    FLUXES.maps.(varnames{p}) = NaN(ncells,1);
    
    for k=1:ncells

        FLUXES.maps.(varnames{p})(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
        % Currently only takes the first soil moisture layer. Modify to
        % split up each soil moisture layer into its own variable.
        
    end
    
end

end