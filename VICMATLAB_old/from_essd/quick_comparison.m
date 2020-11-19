% Quick comparison
%
% Quickly compare NetCDF runoff/baseflow output from VIC image driver with
% USGS gauge data

fluxfile = '/Volumes/HD4/SWOTDA/Data/Tuolumne/WB/fluxes.1999-01-01.nc';
info = ncinfo(fluxfile);
runoff = ncread(fluxfile, 'OUT_RUNOFF');
baseflow = ncread(fluxfile, 'OUT_BASEFLOW');

timevector = datetime(1999,1,1):days(1):datetime(1999,12,31);

mean_runoff = squeeze(nanmean(nanmean(runoff),2));
mean_baseflow = squeeze(nanmean(nanmean(baseflow),2));

combined_flow = mean_runoff + mean_baseflow; % mm/day

% Rough unit conversion
rect_area = 26134; % square kilometers
basin_area = 4850592028; % square meters
basin_area_km2 = 4850592028/1000^2; % Upper Tuolumne is 4850 square km
combined_flow_cms = combined_flow*rect_area*1000/(24*3600);

figure, subplot(2,1,1)
jsplot(timevector, combined_flow, 'Combined runoff and baseflow', 'Time', 'Flow (m^3/s)', 16)

% USGS gauge near the outlet of the basin
gage_file = '/Users/jschap/Documents/Classes/GEOG207 Seminar/Data/QMIN/streamflow_data/archived_data/data/11274790.txt';
dat = dlmread(gage_file, ',', 1, 0);
gagetime = datetime(dat(:,1), dat(:,2), dat(:,3));
gagevals = dat(:,4); % cfs
gagevals_cms = gagevals*(12/39.37)^3;
subplot(2,1,2)
jsplot(gagetime, gagevals, 'USGS gauge 1127490', 'Time', 'Flow (m^3/s)', 16)

