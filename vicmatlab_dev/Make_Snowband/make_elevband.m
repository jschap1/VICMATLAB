% Make elevation band file
%
% Based on elevband.c code from Laura Bowling 8/16/1998
% Written 4/2/2019 JRS
% 
% Could vectorize to speed up, but doesn't seem worth it right now
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

function L = make_elevband(soils, dem_name, coarse_res, numbands, min_delta, min_fract, outfile, verbose)

[dem, R] = geotiffread(dem_name);

if min(dem(:)) < -1000
    nodataval = min(dem(:));
    dem(dem==nodataval) = NaN;
end
% assumes the no data value is some large negative number

% figure, imagesc(dem)

if abs(R.CellExtentInLatitude - R.CellExtentInLongitude) > 1e-4
    error('resolution is different in x and y directions')
end

fine_res = R.CellExtentInLatitude;

% number of fine resolution cells in a coarse resolution cell (approximately)
% numcells = floor(coarse_res/fine_res); 
numcells = ceil(coarse_res/fine_res); 

ncells = size(soils,1); % number of coarse resolution grid cells
cellnum = soils(:,2);
lat = soils(:,3);
lon = soils(:,4);
ymin = R.LatitudeLimits(1);
xmin = R.LongitudeLimits(1);
rows = R.RasterSize(1); % number of rows in the fine resolution DEM
cols = R.RasterSize(2);  

L = zeros(ncells, 1+3*numbands);

% loop over coarse resolution grid cells
for k=1:ncells
    
    % the indexing seems more or less right
    % might be some issues near the edge, though
    % it would be worth plotting the snowband file afterward to check
    % everything is as expected 
    
    corner_lat = lat(k) + coarse_res/2;
    corner_lon = lon(k) - coarse_res/2;
   
%     corner_lat = lat(k);
%     corner_lon = lon(k);
    
    corner_lat_ind = rows - round((corner_lat - ymin)/fine_res);
    corner_lon_ind = round((corner_lon - xmin)/fine_res);
    
    % accounts for rounding error at the edges
    if corner_lon_ind <= 0
        corner_lon_ind = 1;
    end
    
    if corner_lat_ind <= 0
        corner_lat_ind = 1;
    end
    
    if corner_lon_ind > size(dem,2) - numcells
        corner_lon_ind = size(dem,2) - numcells;
    end
    
    if corner_lat_ind > size(dem,1) - numcells
        corner_lat_ind = size(dem,1) - numcells;
    end
           
    % check
    
%     figure, imagesc([-180,180],[-90,90],flipud(dem)), 
%     set(gca, 'ydir', 'normal')
%     hold on
%     plot(corner_lon, corner_lat, 'ro')
% 
% figure, imagesc([-180,180],[-90,90],dem)
% 
%     figure, imagesc(dem), hold on
%     plot(i_corner, j_corner, 'ro')
%     plot(corner_lon_ind+numcells, corner_lat_ind+numcells, 'ro')
        
    % initialize
    minval = 9999;
    maxval = 0;
    sum1 = 0;
    count = 0;
    numzero = 0;
    
%     dem(corner_lat_ind, corner_lon_ind)
        
    % loop over the fine resolution grid cells within the coarse res cell
    for i=corner_lat_ind:corner_lat_ind+numcells
        for j=corner_lon_ind:corner_lon_ind+numcells
            
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
    


%% If there are no fine-resolution data (due to georeferencing error), set one elevation band
    if count == 0
        disp(['Error: no data in cell ', num2str(k)])
        disp(['Lat: ', num2str(lat(k)), '; Lon: ', num2str(lon(k))])
        
        area_fract = zeros(numbands,1);
        area_fract(1) = 1;
        avg_elev = zeros(numbands,1);
        avg_elev(1) = soils(k,5);
        pfactor = area_fract;
        L(k,:) = [cellnum(k), area_fract', avg_elev', pfactor'];
        
%% Else (in the majority of cases), use the fine-resolution elevations to make elevation bands       
    else
        meanval = sum1/count; % mean elevation value within the coarse res cell
        
            % calculate difference in elevation, used to determine partitioning
    % into elevation bands
    if (maxval-meanval) > (meanval-minval)
        diff1 = maxval-meanval;
    else
        diff1 = meanval-minval;
    end
    % diff1 = the smaller of the difference between the mean and the max or
    % the mean and the min elevation values w/in the coarse grid cell
    
    % creating elevation increments
    increment = zeros(numbands, 1);
    inc_count = zeros(numbands, 1);
    inc_sum = zeros(numbands, 1);
    increment(1) = meanval - diff1;
    for m=2:numbands
        increment(m) = increment(m-1) + 2*diff1/numbands; % might want to fine-tune this line
%         inc_count(m-1) = 0;
%         inc_sum(m-1) = 0;
    end
    
    if verbose
        if minval<increment(1)
            disp('range higher than min')
        end
        if maxval > increment(numbands)
            disp('warning: range less than maxval for ')
            disp(['Lat: ', num2str(lat(k)), ', Lon: ', num2str(lon(k))])
        end
    end
    increment(numbands) = increment(numbands) + 1; % shifting by 1 to avoid negative values due to roundoff errors
    
    % allocate elevations to bins
    zvals = zeros(count, 1); % get the elevation values for debugging purposes
    ind = 1;
    for i=corner_lat_ind:corner_lat_ind+numcells
        for j=corner_lon_ind:corner_lon_ind+numcells
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
    if verbose
        if count ~= sum(inc_count)
            disp('Not all elevations allocated to bins');
        end
    end
    
    %% Removing bins that are not needed
    
    % check for bins within min_delta
    num_remove = 0;    
    for m=2:numbands
        if (inc_count(m)>0) && (inc_count(m-1) > 0) && (inc_sum(m)/inc_count(m) - inc_sum(m-1)/inc_count(m-1) < min_delta)
            % combine similar bins
            inc_count(m) = inc_count(m) + inc_count(m-1);
            inc_sum(m) = inc_sum(m) + inc_sum(m-1);
            inc_count(m-1) = 0;
            inc_sum(m-1) = 0;
            num_remove = num_remove + 1;
        end
    end
    
    % Check for bins less than min_fract
    for m=1:numbands-1
        fract = inc_count(m)/count;
        if fract <= min_fract && inc_count(m) > 0 && inc_count(m+1) > 0
            inc_count(m+1) = inc_count(m) + inc_count(m+1);
            inc_sum(m+1) = inc_sum(m) + inc_sum(m+1);
            inc_count(m) = 0;
            inc_sum(m) = 0;
            num_remove = num_remove + 1;
        end
    end
    
    % handle the last elevation band separately (m=numbands)
    if inc_count(numbands)/count < min_fract && inc_count(numbands) > 0 && inc_count(numbands-1) > 0
        inc_count(numbands-1) = inc_count(numbands-1) + inc_count(numbands);
        inc_sum(numbands-1) = inc_sum(numbands-1) + inc_sum(numbands);
        inc_count(numbands) = 0;
        inc_sum(numbands) = 0;
        num_remove = num_remove + 1;
    end
    
    
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
    
    % Make sure nothing is negative
    zero_ind = find(L(k,:) < 0);
    L(k, zero_ind) = 0;
    % Negative snow band elevations cause VIC to throw an error
    
    if mod(k,1000)==0 % counter to display progress
        disp(round(100*k/ncells, 3));
    end    
        
    end
    
end

% try to ensure that numbers are output w/sufficient precision
precis = numel(num2str(max(cellnum)));
dlmwrite(outfile, L, 'Delimiter', '\t','precision', precis+1);

return