function FLUXES = basinavgts(FLUXES)

ncells = length(fieldnames(FLUXES.ts));
cellnames = fieldnames(FLUXES.ts);
varnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;

% triple loop = bad. optimize.
for p = 1:length(varnames)
    
    for t = 1:length(FLUXES.time)
        
        sum1 = 0;
        
        for k=1:ncells
            sum1 = sum1 + FLUXES.ts.(cellnames{k}).(varnames{p})(t);   
        end
        
        FLUXES.avgts.(varnames{p})(t) = sum1/ncells;
        
    end
    
end
    
return