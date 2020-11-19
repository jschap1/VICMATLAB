% THIS HAS BEEN MIGRATED TO R BECAUSE R IS BETTER FOR GIS
%
% Fine resolution basin mask (at the DEM resolution)
% Coarse resolution basin mask (at the VIC model resolution)

[fine,R] = arcgridread('/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.asc');
[coarse,R] = arcgridread('/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah_coarse.asc');
 
% number of fine grid cells within one coarse grid cell (64)
coarse_res = 1/8;
fine_res = 1/64;
nfine = (coarse_res/fine_res)^2;

% Compute fraction file

% Coarse and fine rasters must be such that the number of fine resolution
% grid cells within one coarse resolution grid cell is a perfect square

fract = NaN(size(coarse));
for i=1:11
    for j=1:18
        count = 0;
        for m=1+sqrt(nfine)*(i-1):sqrt(nfine)*i
            for n=1+sqrt(nfine)*(j-1):sqrt(nfine)*j
                count = count + fine(m,n);
            end
        end
        fract(i,j) = count;
    end
end

% The above does not work...


i=1;j=1;
count = 0;
for m=1:sqrt(nfine)
    for n=1:sqrt(nfine)
        count = count + fine(m,n);
    end
end
fract(i,j) = count;
        
% General







inpolygon(basin.Y, basin.Y)

% Name and location of basin shapefile
shpname = '/Users/jschapMac/Desktop/UpperKaweah/upper_kaweah.shp';

% Name and location of soil parameter file
soilname = '/Users/jschapMac/Desktop/UpperKaweah/soils_spec.KAWEAH';

basin = shaperead(shpname); % polygon
soil = load(soilname); % grid cells

slat = soil(:,3);
slon = soil(:,4);

% This makes a raster from X Y data.
[Z,R] = vec2mtx(basin.Y, basin.X, 1);

% Compute the portion of the grid cell that is in the basin:
