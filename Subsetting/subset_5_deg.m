% Subsetting the domain and parameter files into * degree tiles
%
% Filename includes the lower-left coordinates of the tile

addpath('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting/final')
global_domain = '/Volumes/HD_ExFAT/output/VICGlobal_domain.nc';
global_params = '/Volumes/HD_ExFAT/output/VICGlobal_params.nc';

res = 10; % set the desired resolution for subsetting

% Subset domain
for i=1:(180/res)
    for j=1:(360/res)   
        extent{i,j} = horzcat([-180 + res*(j-1); -180 + res*j], [-90 + res*(i-1); -90 + res*i]);
        % Subset domain
        outname_domain = ['/Volumes/HD_ExFAT/output/domain_tiles/domain_' ...
            'lon_' num2str(-180 + 10*(j-1)) '_lat_' num2str(-90 + 10*(i-1)) '.nc'];
        subset_domain(extent{i,j}, global_domain, outname_domain)
    end
end

% Subset parameter
for i=1:(180/res)
    for j=1:(360/res)  
        extent{i,j} = horzcat([-180 + res*(j-1); -180 + res*j], [-90 + res*(i-1); -90 + res*i]);      
        outname_params = ['/Volumes/HD_ExFAT/output/parameter_tiles/params_' ...
            'lon_' num2str(-180 + res*(j-1)) '_lat_' num2str(-90 + res*(i-1)) '.nc'];
        subset_parameter(extent{i,j}, global_params, outname_params)  
    end
end

