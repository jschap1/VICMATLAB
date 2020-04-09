% Make cropping rectangle
%
% OUTPUTS
% rect = cropping rectangle for imcrop
% out_lon = longitude values for the cropping extent (vector)
% out_lat = latitude values for the cropping extent (vector)

function [rect, out_lon, out_lat, height, width] = make_cropping_rectangle(lon, lat, lon_range, lat_range, xres, yres) 

K = 4000; % parameter
M = 50; % another parameter

% must use grid cell centers for makerefmat
R1 = makerefmat(min(lon), min(lat), xres, yres);

%     this formulation guarantees whole number indices
%     [ymin, xmin] = latlon2pix(R1, min(lat)+yres*10, min(lon)+xres*5) % lat, lon

% find appropriate minimum values for the cropping rectangle
minval_opt_1 = M;
minval_opt_2 = M;
for p=1:K
    tmp1 = abs(min(lon)+ p*xres - lon_range(1));
    tmp2 = abs(min(lat)+ p*yres - lat_range(1));
    if tmp1 < minval_opt_1
        minval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < minval_opt_2
        minval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
minlon = min(lon)+ p_opt_1*xres;
minlat = min(lat)+ p_opt_2*yres;

% make a 1-pixel border just to be safe
minlon = minlon - xres;
minlat = minlat - yres;

% find appropriate maximum values for the cropping rectangle
maxval_opt_1 = M;
maxval_opt_2 = M;
for p=1:K
    tmp1 = abs(minlon + p*xres - lon_range(2));
    tmp2 = abs(minlat + p*yres - lat_range(2));
    if tmp1 < maxval_opt_1
        maxval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < maxval_opt_2
        maxval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
maxlon = minlon + p_opt_1*xres;
maxlat = minlat + p_opt_2*yres;  

% make a 1-pixel border just to be safe
maxlon = maxlon + xres;
maxlat = maxlat + yres;

[ymin, xmin] = latlon2pix(R1, minlat, minlon);
[ymax, xmax] = latlon2pix(R1, maxlat, maxlon);
width = xmax - xmin;
height = ymax - ymin;

%     rect = [floor(xmin), floor(ymin), ceil(width), ceil(height)];
rect = [xmin, ymin, width, height];

% I have to do this to get the code to work...
height = height+1;
width = width+1;

out_lat = minlat:yres:maxlat;
out_lon = minlon:xres:maxlon;

return


