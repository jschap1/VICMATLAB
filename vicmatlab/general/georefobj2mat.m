% Georeferencing object to georeferencing matrix
%
% Written for latlon coordinates
%
% INPUTS
% R = georeferencing object
% type = "LL" or "center"

function Rmat = georefobj2mat(R, type)

xres = R.CellExtentInLongitude;
yres = R.CellExtentInLatitude;

switch type
    case 'LL'
        xmin = R.LongitudeLimits(1) + xres/2;
        ymin = R.LatitudeLimits(1) + yres/2;
    case 'center'
        xmin = R.LongitudeLimits(1);
        ymin = R.LatitudeLimits(1);
    otherwise
        disp('Please provide information about raster format')
end

Rmat = makerefmat(xmin, ymin, xres, yres);


return