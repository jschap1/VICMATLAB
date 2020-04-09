function bndpts = GetBndPoints(fdir, river)

% Calculates upstream boundary points using river and flow direction
% rasters as input. Part of the LISFLOOD input file preparation workflow
%
% I would use Kostas' Python function to do this, but it is thrown off by
% my rasters. Perhaps it would work if I created all the inputs in GRASS.
%
% Actually, this function diverges from Kostas (correct) method. It ought
% to check if river pixels, not just any pixels, are flowing into the
% current cell. DO NOT USE AS IS.
%
% Example usage
% cd /Users/jschapMac/Desktop/Tuolumne/TuoSub/GIS
% addpath TempGIS
% [fdir, R] = geotiffread('fdir.tif');
% [river, r] = arcgridread('river.asc');
% bndpts = GetBndPoints(fdir, river)
% csvwrite('bndpts.csv', bndpts)

bndpts = [];
[nrow, ncol] = size(fdir);

for i=1:nrow
    for j=1:ncol
        
        if river(i,j) == 1
            up_cells = CheckCells(fdir, i, j);
            if isempty(up_cells)
                bndpts = [bndpts; [i,j]];
            end
        end
        
    end
end

% Plotting to check results
figure, imagesc(river), hold on
plot(bndpts(:,2),bndpts(:,1),'r*')

fd = double(fdir);
fd(fd<0) = NaN;
figure, imagesc(fd), hold on
plot(bndpts(:,2),bndpts(:,1),'r*')

return

function up_cells = CheckCells(fdir, i, j)

    % Checks adjacent cells and returns a list of cell indices for any
    % cells draining into cell(i,j)
    
    up_cells = [];
    
    try % deals with edge cells
    if fdir(i,j+1) == 7 % check if the cell to the east is flowing west
        up_cells = [up_cells; [i, j+1]];
    end
    catch
    end
    
    try
    if fdir(i,j-1) == 3 % west to east
        up_cells = [up_cells; [i, j-1]];
    end
    catch
    end
    
    try
    if fdir(i+1,j) == 1 % south to north
        up_cells = [up_cells; [i+1, j]];
    end
    catch
    end
    
    try
    if fdir(i-1,j) == 5 % north to south
        up_cells = [up_cells; [i-1, j]];
    end
    catch 
    end
    
    try
    if fdir(i+1,j+1) == 8 %
        up_cells = [up_cells; [i+1, j+1]];
    end
    catch
    end
    
    try
    if fdir(i+1,j-1) == 2 % 
        up_cells = [up_cells; [i+1, j-1]];
    end
    catch
    end
    
    try
    if fdir(i-1,j+1) == 6 % 
        up_cells = [up_cells; [i-1, j+1]];
    end
    catch
    end
    
    try
    if fdir(i-1,j-1) == 4 % 
        up_cells = [up_cells; [i-1, j-1]];
    end
    catch
    end
    
return