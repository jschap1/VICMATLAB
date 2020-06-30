% Geotiffread2 
%
% Modified version of geotiffread that returns lat and lon coordinates

function [A, R, lon, lat] = geotiffread2(fname)

try
    [A, ~, R] = geotiffread(fname);
    A = flipud(A);
    Rmat = georefobj2mat(R, 'LL');
catch
    [A, R] = geotiffread(fname);
    A = flipud(A);
    Rmat = georefobj2mat(R, 'LL');    
end

[lon, lat] = pixcenters(Rmat, size(A));

% Remove no data values
A(A < -1000) = NaN;

return