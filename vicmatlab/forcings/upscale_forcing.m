% Upscale forcing
%
% Upscales forcing data from one resolution to another

function upscale_forcing(forcdir, prefix, newres, oldres)

[lon, lat] = get_coordinates_from_VIC_file(forcdir, prefix);



lat = newsoils(:,3);
lon = newsoils(:,4);
f = newres/oldres;

% Upscale similar to upscale_vegpars or upscale_soils

runcell = xyz2grid(lon, lat, oldsoils(:,1));
[m1, n1] = size(runcell); % original dimensions

% new dimensions
m2 = m1/f;
n2 = floor(n1/f);

% Upscale (or downscale) (method from Walter Roberson in MATLAB Answers)
minlat = min(lat);
maxlat = max(lat);
minlon = min(lon);
maxlon = max(lon);
oldlats = linspace(minlat, maxlat, m1);
oldlons = linspace(minlon, maxlon, n1)';
newlats = linspace(minlat, maxlat, m2);
newlons = linspace(minlon, maxlon, n2)'; % resolution is a couple hundred meters off....

runcell2 = interp2(oldlats, oldlons, runcell', newlats, newlons, 'linear');
runcell2 = runcell2';
runcell2(isnan(runcell2))=0;
runcell3 = runcell2(:); % indexing

return