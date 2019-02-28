function FLUXES = timeavgfluxmaps(FLUXES)

ncells = length(fieldnames(FLUXES.ts));
cellnames = fieldnames(FLUXES.ts);
varnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;

for p = 1:length(varnames)
    
    FLUXES.avgmaps.(varnames{p}) = NaN(ncells,1);
    if strcmp(varnames{p}, 'moist')
        nlayers = size(FLUXES.ts.(cellnames{1}).(varnames{p}),2);
        FLUXES.avgmaps.(varnames{p}) = NaN(ncells,nlayers);
    end
    
    for k=1:ncells
        FLUXES.avgmaps.(varnames{p})(k,:) = mean(FLUXES.ts.(cellnames{k}).(varnames{p}));
    end
    
end

return