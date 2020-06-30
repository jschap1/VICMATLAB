% Get coordinates from VIC file
%
% INPUTS
% folder = name of directory where VIC files are stored
% option = file prefix ('Forcing' or 'Flux', for example)
%
% OUTPUTS
% lat lon coordinates for each VIC file

function [lon, lat] = get_coordinates_from_VIC_file(folder, prefix)

fnames = dir(fullfile(folder, [prefix '*']));
ncells = length(fnames);
lat = zeros(ncells,1);
lon = zeros(ncells,1);

n_underscores = length(strfind(prefix, '_'));

for k=1:ncells

    % one or the other code block is used, depending on the format of the
    % filenames
    if n_underscores == 1
      tmp = strsplit(fnames(k).name, '_');
      lat(k) = str2double(tmp{2});
      tmp = strsplit(tmp{3}, '.txt');
      lon(k) = str2double(tmp{1});
    elseif n_underscores == 2
      tmp = strsplit(fnames(k).name, '_');
      lat(k) = str2double(tmp{3});
      tmp = strsplit(tmp{4}, '.txt');
      lon(k) = str2double(tmp{1});
    end
end

return
