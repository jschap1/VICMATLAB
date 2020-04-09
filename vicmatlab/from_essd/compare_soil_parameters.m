% Plots soil parameters for UCRB and UMRB
% Compare Livneh and VICGlobal parameters
%
% Written 1/27/2020 JRS
%
% Sample inputs:
%
% % Upper Mississippi
% L15_dir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/L15/soils'
% vg_dir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/soils'
% 
% % Upper Colorado
% vg_dir = '/Volumes/HD4/SWOTDA/Data/Colorado/soils';
% L15_dir = '/Volumes/HD4/SWOTDA/Data/Colorado/L15/soils';
% outdir = '/Users/jschap/Documents/Research/VICGlobal/Figures/';

function compare_soil_parameters(vg_dir, L15_dir, outdir)

bounds = struct();

% UMRB
% bounds.precip = [400,1200];
% bounds.temp = [0,15];
% bounds.d1 = [0, 0.5];
% bounds.d2 = [0, 2];
% bounds.d3 = [0, 2];
% bounds.elev = [0, 600];
% bounds.b = [0, 0.5];
% bounds.ds = [0, 1];
% bounds.ws = [0, 1];
% bounds.dsmax = [0, 40];

% UCRB
bounds.precip = [0,1000];
bounds.temp = [-5,15];
bounds.d1 = [0, 0.5];
bounds.d2 = [0, 2];
bounds.d3 = [0, 2];
bounds.elev = [0, 4000];
bounds.b = [0, 0.2];
bounds.ds = [0, 0.4];
bounds.ws = [0.2, 1];
bounds.dsmax = [0, 120];

vicglobal = struct();
livneh = struct();

h1 = 1850;
w1 = 850;

%% Figure A 

[vicglobal.annP, R, lon, lat] = geotiffread2(fullfile(vg_dir, 'annual_prec.tif'));
vicglobal.avgT = geotiffread2(fullfile(vg_dir, 'avg_T.tif'));
vicglobal.d1 = geotiffread2(fullfile(vg_dir, 'depth1.tif'));
vicglobal.d2 = geotiffread2(fullfile(vg_dir, 'depth2.tif'));
vicglobal.d3 = geotiffread2(fullfile(vg_dir, 'depth3.tif'));

livneh.annP = geotiffread2(fullfile(L15_dir, 'annual_prec.tif'));
livneh.avgT = geotiffread2(fullfile(L15_dir, 'avg_T.tif'));
livneh.d1 = geotiffread2(fullfile(L15_dir, 'depth1.tif'));
livneh.d2 = geotiffread2(fullfile(L15_dir, 'depth2.tif'));
livneh.d3 = geotiffread2(fullfile(L15_dir, 'depth3.tif'));

figure

set(gcf, 'Position',  [100, 100, h1, w1])

subplot(2,5,1)
plotraster(lon, lat, vicglobal.annP, 'Annual prec. (mm)', 'Lon', 'Lat')
caxis(bounds.precip)

subplot(2,5,2)
plotraster(lon, lat, vicglobal.avgT, 'Average temp. (deg. C)', 'Lon', 'Lat')
caxis(bounds.temp)

subplot(2,5,3)
plotraster(lon, lat, vicglobal.d1, 'L1 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d1)

subplot(2,5,4)
plotraster(lon, lat, vicglobal.d2, 'L2 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d2)

subplot(2,5,5)
plotraster(lon, lat, vicglobal.d3, 'L3 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d3)

subplot(2,5,6)
plotraster(lon, lat, livneh.annP, 'Annual prec. (mm)', 'Lon', 'Lat')
caxis(bounds.precip)

subplot(2,5,7)
plotraster(lon, lat, livneh.avgT, 'Average temp. (deg. C)', 'Lon', 'Lat')
caxis(bounds.temp)

subplot(2,5,8)
plotraster(lon, lat, livneh.d1, 'L1 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d1)

subplot(2,5,9)
plotraster(lon, lat, livneh.d2, 'L2 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d2)

subplot(2,5,10)
plotraster(lon, lat, livneh.d3, 'L3 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d3)

saveas(gcf, fullfile(outdir, 'soil_par_comparison_plot_A.png'));

%% Figure B

vicglobal.elev = geotiffread2(fullfile(vg_dir, 'elev.tif'));
vicglobal.b = geotiffread2(fullfile(vg_dir, 'b_infilt.tif'));
vicglobal.ds = geotiffread2(fullfile(vg_dir, 'ds.tif'));
vicglobal.ws = geotiffread2(fullfile(vg_dir, 'ws.tif'));
vicglobal.dsmax = geotiffread2(fullfile(vg_dir, 'dsmax.tif'));

livneh.elev = geotiffread2(fullfile(L15_dir, 'elev.tif'));
livneh.b = geotiffread2(fullfile(L15_dir, 'infilt.tif'));
livneh.ds = geotiffread2(fullfile(L15_dir, 'ds.tif'));
livneh.ws = geotiffread2(fullfile(L15_dir, 'ws.tif'));
livneh.dsmax = geotiffread2(fullfile(L15_dir, 'dsmax.tif'));

figure

set(gcf, 'Position',  [100, 100, h1, w1])

subplot(2,5,1)
plotraster(lon, lat, vicglobal.elev, 'Elevation (m)', 'Lon', 'Lat')
caxis(bounds.elev)
subplot(2,5,2)
plotraster(lon, lat, vicglobal.b, 'VIC parameter (b)', 'Lon', 'Lat')
caxis(bounds.b)
subplot(2,5,3)
plotraster(lon, lat, vicglobal.ds, 'Ds', 'Lon', 'Lat')
caxis(bounds.ds)
subplot(2,5,4)
plotraster(lon, lat, vicglobal.ws, 'Ws', 'Lon', 'Lat')
caxis(bounds.ws)
subplot(2,5,5)
plotraster(lon, lat, vicglobal.dsmax, 'Dsmax (mm/day)', 'Lon', 'Lat')
caxis(bounds.dsmax)

subplot(2,5,6)
plotraster(lon, lat, livneh.elev, 'Elevation (m)', 'Lon', 'Lat')
caxis(bounds.elev)
subplot(2,5,7)
plotraster(lon, lat, livneh.b, 'VIC parameter (b)', 'Lon', 'Lat')
caxis(bounds.b)
subplot(2,5,8)
plotraster(lon, lat, livneh.ds, 'Ds', 'Lon', 'Lat')
caxis(bounds.ds)
subplot(2,5,9)
plotraster(lon, lat, livneh.ws, 'Ws', 'Lon', 'Lat')
caxis(bounds.ws)
subplot(2,5,10)
plotraster(lon, lat, livneh.dsmax, 'Dsmax (mm/day)', 'Lon', 'Lat')
caxis(bounds.dsmax)

saveas(gcf, fullfile(outdir, 'soil_par_comparison_plot_B.png'));

%% Figure C (for ESSD paper)

figure

set(gcf, 'Position',  [100, 100, h1, w1])

subplot(2,5,1)
plotraster(lon, lat, vicglobal.b, 'VIC parameter (b)', 'Lon', 'Lat')
caxis(bounds.b)
subplot(2,5,2)
plotraster(lon, lat, vicglobal.ds, 'Ds', 'Lon', 'Lat')
caxis(bounds.ds)
subplot(2,5,3)
plotraster(lon, lat, vicglobal.ws, 'Ws', 'Lon', 'Lat')
caxis(bounds.ws)
subplot(2,5,4)
plotraster(lon, lat, vicglobal.dsmax, 'Dsmax (mm/day)', 'Lon', 'Lat')
caxis(bounds.dsmax)
subplot(2,5,5)
plotraster(lon, lat, vicglobal.d3, 'L3 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d3)

subplot(2,5,6)
plotraster(lon, lat, livneh.b, 'VIC parameter (b)', 'Lon', 'Lat')
caxis(bounds.b)
subplot(2,5,7)
plotraster(lon, lat, livneh.ds, 'Ds', 'Lon', 'Lat')
caxis(bounds.ds)
subplot(2,5,8)
plotraster(lon, lat, livneh.ws, 'Ws', 'Lon', 'Lat')
caxis(bounds.ws)
subplot(2,5,9)
plotraster(lon, lat, livneh.dsmax, 'Dsmax (mm/day)', 'Lon', 'Lat')
caxis(bounds.dsmax)
subplot(2,5,10)
plotraster(lon, lat, livneh.d3, 'L3 soil thickness (mm)', 'Lon', 'Lat')
caxis(bounds.d3)

saveas(gcf, fullfile(outdir, 'soil_par_comparison_plot_C.png'));

return
