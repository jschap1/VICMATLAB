% Get coordinates (lat/lon vectors) from R and dem
function [ulon, ulat] = getlatlon(A, R)

    [nx, ny] = size(A);
    
    
    x=(1:nx)';
    y=ones(length(x),1);
    
    xx = repmat(x, ny, 1);
   
    % For this bit of code, credit to:
    % https://www.mathworks.com/matlabcentral/answers/16270-create-vector-of-repeating-elements-sort-of
    v=repmat(1:ny, [nx 1]);
    yy = v(:); 
    
    [lat, lon] = intrinsicToGeographic(R,xx,yy);

    ulat = unique(lat);
    ulon = unique(lon);
    
end