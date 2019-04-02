% Make elevation band file
%
% Based on elevband.c code from Laura Bowling 8/16/1998
% Written 4/2/2019 JRS
% 
% There are some bugs to work out.
% Could vectorize to speed up
%
% INPUTS
% Soil parameter file (specifically, the first four columns)
% Fine resolution DEM
% numbands = 5 is a good value
%     bands less than min_delta meters are are merged together
%     bands covering less than min_fract of the grid cell area are removed
%
% OUTPUTS
% Elevation band file

function L = make_elevband(soils, dem_name, coarse_res, numbands, min_delta, min_fract, outfile)

[dem, R] = geotiffread(dem_name);

if min(dem(:)) < -1000
    nodataval = min(dem(:));
    dem(dem==nodataval) = NaN;
end
% assumes the no data value is some large negative number

if abs(R.CellExtentInLatitude - R.CellExtentInLongitude) > 1e-4
    error('resolution is different in x and y directions')
end

fine_res = R.CellExtentInLatitude;

% number of fine resolution cells in a coarse resolution cell (approximately)
numcells = floor(coarse_res/fine_res); 

ncells = size(soils,1); % number of coarse resolution grid cells
cellnum = soils(:,2);
lat = soils(:,3);
lon = soils(:,4);
ymin = R.LatitudeLimits(1);
xmin = R.LongitudeLimits(1);
rows = R.RasterSize(1); % number of rows in the fine resolution DEM
    
L = zeros(ncells, 1+3*numbands);

% loop over coarse resolution grid cells
for k=1:ncells
    
    corner_lat = lat(k) + coarse_res/2;
    corner_lon = lon(k) - coarse_res/2;
    
    i_corner = round(rows - ((corner_lat - ymin)/fine_res) - 1);
    j_corner = round((corner_lon - xmin)/fine_res);
    
    % check
%     figure, imagesc(dem), hold on
%     plot(i_corner, j_corner, 'ro')
    
    % initialize
    minval = 9999;
    maxval = 0;
    sum1 = 0;
    count = 0;
    numzero = 0;
    
    % loop over the fine resolution grid cells within the coarse res cell
    for i=i_corner:i_corner+numcells
        for j=j_corner:j_corner+numcells
            
            if ~isnan(dem(i,j))
                
                sum1 = sum1 + dem(i,j);
                count = count + 1;
                
                if dem(i,j) < minval
                    minval = dem(i,j);
                end
                if dem(i,j) > maxval
                    maxval = dem(i,j);
                end
            end
        end
    end
    
%     count = (numcells+1)^2; % should be equivalent
    
    if count == 0
        disp(['Error: no data in cell', num2str(k)])
    else
        meanval = sum1/count; % mean elevation value within the coarse res cell
    end
    
    % calculate difference in elevation, used to determine partitioning
    % into elevation bands
    if (meanval-meanval) > (meanval-minval)
        diff1 = maxval-meanval;
    else
        diff1 = meanval-minval;
    end
    % diff1 = the smaller of the difference between the mean and the max or
    % the mean and the min elevation values w/in the coarse grid cell
    
    % initializing increments (not sure if this is correct)
    increment = zeros(numbands, 1);
    inc_count = zeros(numbands, 1);
    inc_sum = zeros(numbands, 1);
    
    increment(1) = meanval - diff1;
    if minval<increment(1)
        disp('range higher than min')
    end
    
    for m=2:numbands
        increment(m) = increment(m-1) + 2*diff1/numbands;
%         inc_count(m-1) = 0;
%         inc_sum(m-1) = 0;
    end
    
    if maxval > increment(numbands)
        disp('range less than maxval')
    end
    increment(numbands) = increment(numbands) + 1; % shifting by 1 to avoid negative values due to roundoff errors
    
    % allocate elevations to bins
    zvals = zeros(count, 1); % get the elevation values for debugging purposes
    ind = 1;
    for i=i_corner:i_corner+numcells
        for j=j_corner:j_corner+numcells
            if ~isnan(dem(i,j))    
                 
                m = 1;
                if dem(i,j) >= 0 && dem(i,j) < increment(m+1)
                    inc_count(m) = inc_count(m) + 1;
                    inc_sum(m) = inc_sum(m) + dem(i,j);
                end
                
                for m=2:numbands-1
                    if (dem(i,j) >= increment(m)) && (dem(i,j) < increment(m+1))
                        inc_count(m) = inc_count(m) + 1;
                        inc_sum(m) = inc_sum(m) + dem(i,j);
                    end
                end
                
                m = numbands;
                if dem(i,j) >= increment(m)
                    inc_count(m) = inc_count(m) + 1;
                    inc_sum(m) = inc_sum(m) + dem(i,j);
                end
                
                zvals(ind) = dem(i,j);
                ind = ind + 1;
                
            end
        end
    end
    
    % check that all elevations were binned
    if count ~= sum(inc_count)
        disp('Not all elevations allocated to bins');
    end
    
    
    %% Removing bins that are not needed
    
%     % check for bins within min_delta
%     num_remove = 0;
%     for m=2:numbands
%         if (inc_count(m)>0) && (inc_count(m-1) > 0) && (inc_sum(m)/inc_count(m) - inc_sum(m-1)/inc_count(m-1) < min_delta)
%             % combine similar bins
%             inc_count(m) = inc_count(m) + inc_count(m-1);
%             inc_sum(m) = inc_sum(m) + inc_sum(m-1);
%             inc_count(m-1) = 0;
%             inc_sum(m-1) = 0;
%             num_remove = num_remove + 1;
%         end
%     end
%     
%     % Check for bins less than min_fract
%     for m=1:numbands
%         fract = inc_count(m)/count;
%         if fract <= min_fract && inc_count(m) > 0 && inc_count(m+1) > 0
%             inc_count(m+1) = inc_count(m) + inc_count(m+1);
%             inc_sum(m+1) = inc_sum(m) + inc_sum(m+1);
%             inc_count(m) = 0;
%             inc_sum(m) = 0;
%             num_remove = num_remove + 1;
%         end
%     end
%     
%     if inc_count(numbands)/count < min_fract && inc_count(numbands) > 0 && inc_count(numbands-1) > 0
%         inc_count(numbands-1) = inc_count(numbands-1) + inc_count(numbands);
%         inc_sum(numbands-1) = inc_sum(numbands-1) + inc_sum(numbands);
%         inc_count(numbands) = 0;
%         inc_sum(numbands) = 0;
%         num_remove = num_remove + 1;
%     end
    
    
    %% Writing output
    
    % fractional area of each band
    area_fract = zeros(numbands, 1);
    for m=1:numbands
        area_fract(m) = inc_count(m)/count;
    end
    
    % mean elevation for each band
    avg_elev = zeros(numbands, 1);
    pixel_mean = zeros(numbands, 1);
    for m=1:numbands
        if inc_count(m)>0
            avg_elev(m) = inc_sum(m)/inc_count(m);
        else
            avg_elev(m) = 0;
            numzero = numzero + 1;
        end
        pixel_mean(m) = avg_elev(m)*area_fract(m);
    end
    
    % fractional precipitation for each band
%     precip_fract = zeros(numbands, 1);
%     avg_precip = zeros(numbands, 1);
%     pfactor = zeros(numbands, 1);
%     for m=1:numbands
%         if inc_count(m)>0
%             avg_precip(m) = inc_sum(m)/inc_count(m);
%         else
%             avg_precip(m) = 1;
%         end
%         precip_fract(m) = inc_count(m)/count;
%         if precip_fract(m)<0 || pixel_mean(m)==0
%             pfactor(m) = 0;
%         else
%             pfactor(m) = precip_fract(m)*avg_precip(m)/pixel_mean(m);
%         end
%     end

    % assuming precipitation is evenly distributed across elevations
    pfactor = area_fract;

    L(k,:) = [cellnum(k), area_fract', avg_elev', pfactor'];
    
end

dlmwrite(outfile, L, '\t');

return