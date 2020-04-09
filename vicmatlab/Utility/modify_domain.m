% Modify domain for VIC image mode
%
% Goal is to get the domain to be consistent across the different input
% files

domainfile = '/Volumes/HD4/SWOTDA/Data/UMRB/umrb_domain.nc';

mask1 = ncread(domainfile, 'mask');
lon1 = ncread('/Volumes/HD4/SWOTDA/Data/UMRB/umrb_domain.nc', 'lon');
lat1 = ncread('/Volumes/HD4/SWOTDA/Data/UMRB/umrb_domain.nc', 'lat');

forc_livneh = ncread('/Volumes/HD4/SWOTDA/Data/UMRB/livneh_forc.1993.nc', 'tas');
forc_livneh = forc_livneh(:,:,1)';

forc_merra2 = ncread('/Volumes/HD4/SWOTDA/Data/UMRB/merra2_forc.1993.nc', 'tas');
forc_merra2 = forc_merra2(:,:,1)';

figure,
plotraster(lon1, lat1, forc_livneh, 'Livneh Temperature', '', '')

figure,
plotraster(lon1, lat1, forc_merra2, 'MERRA-2 Temperature', '', '')

nan_ind = isnan(forc_livneh(:));
mask_livneh = ones(size(forc_livneh));
mask_livneh(nan_ind) = 0;

ncwrite(domainfile, 'mask', mask_livneh');

figure,
plotraster(lon1, lat1, mask_livneh, 'Mask (L15)', '', '')

paramfile = '/Volumes/HD4/SWOTDA/Data/UMRB/umrb_params.nc';
ncwrite(paramfile, 'run_cell', mask_livneh')

run_cell = ncread(paramfile, 'run_cell');