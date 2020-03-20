% Convert .hdf files to .mat files
% Use this function to simplify importing the MERRA-2 data with
% fileDatastore in prepare_merra2_forcings
% merra_dir = '/Volumes/HD_ExFAT/Sample';

function convert_hdf_to_mat(merra_dir)

fnames = dir(fullfile(merra_dir, '*.hdf'));
nfiles = length(fnames);

for k=1:nfiles
    dat = hdfread(fullfile(merra_dir, fnames(k).name), 'EOSGRID', 'Fields', 'T2M');
    save(fullfile(merra_dir, [fnames(k).name, '.mat']))
end

return