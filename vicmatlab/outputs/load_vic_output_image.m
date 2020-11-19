% Loads outputs from VIC image mode into MATLAB
%
% Written 8/4/2020 JRS

function FLUXES = load_vic_output_image(vic_out_name)

% Load NetCDF file contents into Matlab
info = ncinfo(vic_out_name);
numvars = length(info.Variables);
varnames = cell(numvars,1);

for p=1:numvars
    varnames{p} = info.Variables(p).Name;
    var = ncread(vic_out_name, varnames{p});
    
    if ndims(var) == 4
        disp([varnames{p} 'has four dims'])
        FLUXES.(varnames{p}) = permute(var, [2,1,3,4]);
    else
        FLUXES.(varnames{p}) = permute(var, [2,1,3]);
    end
end

% Convert to datetime:
FLUXES.time = datestr(FLUXES.time+367);
FLUXES.time = datetime(FLUXES.time);

return