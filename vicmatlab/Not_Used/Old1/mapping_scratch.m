% Sample map

lat = xyz(:,1);
lon = xyz(:,2);
value = xyz(:,3);

[X,Y] = meshgrid(lat,lon);
Z = NaN(size(X));



tab = table(lat, lon, value);
tab = sortrows(tab);

value2 = reshape(value, length(x), length(y));

imagesc(x,y,value2)
colorbar

%%
lat = [41, 39, 41, 39]';
lon = -[119,119,118,118]';
coords = [lat, lon];
value = [4,2,3,2]';
[X,Y] = meshgrid(lat,lon);
Z = NaN(size(X));

find(coords == [X,Y])

% Replace the NaN values with values from VIC results if the mask is equal
% to one at that cell.

mask = [0,0,0,0;1,1,0,0;1,1,0,0;0,0,0,0];

Z(mask == 1) = value;

