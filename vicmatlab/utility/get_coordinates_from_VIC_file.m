% Get coordinates from VIC file
%
% INPUTS
% folder = name of directory where VIC files are stored
% option = file prefix ('Forcing' or 'Flux', for example)
%
% OUTPUTS
% lat lon coordinates for each VIC file

function [lon, lat] = get_coordinates_from_VIC_file(folder, prefix, method)

fnames = dir(fullfile(folder, [prefix '*']));
ncells = length(fnames);
lat = zeros(ncells,1);
lon = zeros(ncells,1);
for k=1:ncells

    % one or the other code block is used, depending on the format of the
    % filenames
    if strcmp(method, 'forcings')
      tmp = strsplit(fnames(k).name, '_');
      lat(k) = str2double(tmp{2});
      tmp = strsplit(tmp{3}, '.txt');
      lon(k) = str2double(tmp{1});
    elseif strcmp(method, 'fluxes')
      tmp = strsplit(fnames(k).name, '_');
      lat(k) = str2double(tmp{3});
      tmp = strsplit(tmp{4}, '.txt');
      lon(k) = str2double(tmp{1});
    end
end

return
