% Finds row and column indices of a gage for the routing model station 
% location file. Could probably make this more efficient/eliminate 
% looping through all the grid cells

% Written by Jacob Schaperow, Aug. 9, 2017

% Find station location in fdir file
[fdir, R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne/RoutingInputs/fdir_in.asc');

gage = [-121.151    37.6];

min_diff = 1/16; % resolution
% the difference between the gage location and the grid cell coordinates
% will never be greater than the grid resolution

% Load basin mask (should be the same grid as the VIC modeling domain)
%[vicdomain, R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne/RoutingInputs/basinmask_coarse.asc');

% check each cell in the fdir file and see if it is closest tot he stations
[nrows, ncols] = size(fdir);

breakflag = 0; % needed bc it is a pain for Matlab to break out of nested loops
for i=1:nrows
    for j=1:ncols
        [lat, lon] = pix2latlon(R,i,j);
        xdif = abs(lon - gage(:,1));
        ydif = abs(lat - gage(:,2));
        
        if xdif < min_diff && ydif < min_diff
            disp('found it')
            % stnloc file uses col, row, starting at the bottom left of the raster.
            row = nrows - i;
            col = j;
            breakflag = 1;
            break;
        end
        
    end
    
    if breakflag
        break;
    end
    
end

% [col, row] = 14, 10 or 14, 22;