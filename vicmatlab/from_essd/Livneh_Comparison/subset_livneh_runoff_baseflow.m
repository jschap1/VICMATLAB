% Subsets the Livneh et al. (2015) flux outputs to a particular time and
% place
%
% Sums runoff and baseflow together
%
% Written 1/14/2020 JRS
% 
% Modified 1/29/2020 JRS
% Now is a function, not a script

% Inputs
% addpath('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting/Subsetting')
% start_date = datetime(1979, 10, 1); % must start on the first day of the month
% end_date = datetime(2011, 9, 30);
% basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% fluxdir = '/Volumes/HD3/Livneh_2015/Fluxes';
% outdir = '/Volumes/HD4/SWOTDA/Data/Colorado/L15';

function subset_livneh_runoff_baseflow(start_date, end_date, basin_mask, fluxdir, outdir)

%% Do not modify below here

fluxnames = dir(fullfile(fluxdir, '*.nc'));

[mask1, R1, lon1, lat1] = geotiffread2(basin_mask);
figure, plotraster(lon1, lat1, mask1, '','','')

% find the file containing the start date and read it
% continue reading files until you find the end date
% append data as you go

% Initialize arrays
ndays = days(end_date - start_date);
ncells = sum(sum(mask1>0));
runoff_sub = NaN(ndays, ncells);
baseflow_sub = NaN(ndays, ncells);

current_date = start_date;
first_iter = 1; % flag for first iteration
nt_old = 0;
d1 = 1;

while current_date < end_date
    
    current_year = year(current_date);
    current_month = month(current_date);
        
    if current_month < 10
        month_str = ['0' num2str(current_month)]; 
    else
        month_str = num2str(current_month);
    end
    filename = ['Fluxes_Livneh_NAmerExt_15Oct2014.' num2str(current_year) month_str '.nc'];
    filename = fullfile(fluxdir, filename);
    
    if first_iter
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        time_number = ncread(filename, 'time');
    end
    
    runoff = ncread(filename, 'Runoff');
    baseflow = ncread(filename, 'Baseflow');
    
    % Keep only the portion of the output that overlaps the study domain
    % Same methodology as vicinputworkflow for the subsetting    
        
    % Indices for days
    nt = size(runoff, 3);
    d1 = d1 + nt_old;
    d2 = d1 + nt - 1;
    nt_old = nt;

    [runoff_sub(d1:d2,:), tuo_lon, tuo_lat] = subset_netcdf_w_geotiffmask(runoff, lon, lat, basin_mask);
    [baseflow_sub(d1:d2,:), ~, ~] = subset_netcdf_w_geotiffmask(baseflow, lon, lat, basin_mask);
    
    current_date = current_date + nt;
    
    disp(current_date)
    
end

%% Compute the "discharge" time series and plot it

timevector = start_date:end_date;
summed_discharge = nansum(runoff_sub,2) + nansum(baseflow_sub,2);
figure
jsplot(timevector, summed_discharge, 'Summed runoff and baseflow (UCRB)','Time','Q (mm)', 18)

% Convert to m^3/s
A = 35.5; % average area of each grid cell in the Upper Miss. basin
summed_discharge_cms = summed_discharge*A*1000/(24*3600);

figure
jsplot(timevector, summed_discharge_cms, 'Discharge (UMRB)','Time','Q (m^3/s)', 18)

% Report DOWY average values

% outdir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/out_wy1993/processed';
save(fullfile(outdir, 'livneh_q.mat'), 'summed_discharge', 'runoff_sub', 'baseflow_sub', 'timevector');
disp(['Saved ' fullfile(outdir, 'livneh_q.mat')])

return