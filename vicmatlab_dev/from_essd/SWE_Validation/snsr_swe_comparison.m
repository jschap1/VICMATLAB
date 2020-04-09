% Compares SWE between Tuolumne simulation and SNSR reanalysis data

addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB'))

%% Make mask of SNSR SWE

[snsr1, R] = geotiffread('/Volumes/HD3/SNSR/cropped_resampled_tuo/SN_SWE_WY1999_dowy_6_cr.tif');
snsr1 = flipud(snsr1);
snsr1(snsr1<0) = NaN;
Rmat = georefobj2mat(R, 'LL');
[lon, lat] = pixcenters(Rmat, size(snsr1));
figure, plotraster(lon, lat, snsr1, 'SNSR SWE (DOWY 6)', 'Lon', 'Lat')

snsr_mask = snsr1;
snsr_mask(snsr1>=0) = 1;
snsr_mask(isnan(snsr1)) = 0;
figure, plotraster(lon, lat, snsr_mask, 'mask', 'Lon', 'Lat')
% geotiffwrite('/Volumes/HD3/SNSR/snsr_mask2.tif', snsr_mask, Rmat)

%% Calculate VIC-simulated time series (VICGlobal parameters)

wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/Classic_VICGlobal/Raw/wb';
eb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/Classic_VICGlobal/Raw/eb';
timestep_out = 'daily';
info = get_vic_run_metadata(wb_out_dir, eb_out_dir, timestep_out);
swe_vic = read_SWE(info.wb_out_dir, info.ncells, info.nt);

% swe_vic_masked = double(swe_vic.*flipud(snsr_mask_vect));
% swe_vic_map = xyz2grid(flipud(info.lon), info.lat, swe_vic_masked(:,365));
% figure
% plotraster(info.lon, info.lat, swe_vic_map, 'VIC SWE (DOY 365)','Lon','Lat')
% 
% swe_vic2 = (swe_vic);
% % swe_vic2 = flipud(swe_vic);
% lon2 = flipud(info.lon);
% lat2= info.lat;
% swe_vic_map = xyz2grid(lon2, lat2, swe_vic2(:,365));
% figure
% plotraster(info.lon, info.lat, swe_vic_map, 'VIC SWE (DOY 365)','Lon','Lat')
% 
% snsr_mask(snsr_mask==0) = NaN;
% snsr_mask_t = flipud(snsr_mask);
% snsr_mask_vect = snsr_mask_t(:);
% 
% % snsr_mask_t = snsr_mask';
% % snsr_mask_vect = snsr_mask_t(:);

% swe_vic_masked = double(swe_vic.*snsr_mask_vect);

swe_vic_map = fliplr(xyz2grid(info.lon, info.lat, swe_vic(:,365)));
figure
plotraster(info.lon, info.lat, swe_vic_map, 'VIC SWE (DOY 365)','Lon','Lat')

nt = size(swe_vic, 2);
temp_outdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/swe_vg';
mkdir(temp_outdir);
for t=1:nt
    swe_vic_map = fliplr(xyz2grid(info.lon, info.lat, swe_vic(:,t)));
    geotiffwrite(fullfile(temp_outdir, ['swe_vic_' num2str(t) '.tif']), swe_vic_map, R);
end

swe_vic_map_avg = nanmean(swe_vic, 2);
swe_vic_map_avg = fliplr(xyz2grid(info.lon, info.lat, swe_vic_map_avg));
figure
plotraster(info.lon, info.lat, swe_vic_map_avg, 'VIC SWE (average)','Lon','Lat')

vic_dates = timebuilder(1980,1,1,2011,12,31,24);

% writetable(table(vic_dates), fullfile(saveloc, 'vic_dates.txt'))

vic_timevector = datetime(vic_dates(:,1), vic_dates(:,2), vic_dates(:,3));

% Calculate average SWE over SNSR mask
A = swe_vic(snsr_mask(:)==1,:);
swe_vic_avg_ts = nanmean(A,1)';

R_VIC = makerefmat(min(info.lon), min(info.lat), 1/16, 1/16);
swe_mask2 = ones(size(swe_vic_map));
geotiffwrite('/Volumes/HD3/SNSR/swe_mask3.tif', swe_mask2, R_VIC)

% swe_mask2 is used to crop the SNSR SWE in the R script, snsr_swe_conversion2

%% Calculate VIC-simulated time series (Livneh parameters)

wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Raw_EB_FS_1980-2011/wb';
eb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Raw_EB_FS_1980-2011/eb';
timestep_out = 'daily';
info_L15 = get_vic_run_metadata(wb_out_dir, eb_out_dir, timestep_out);
swe_vic_L15 = read_SWE(info_L15.wb_out_dir, info_L15.ncells, info_L15.nt);

% snsr_mask_2 = snsr_mask';
% snsr_mask(isnan(snsr_mask)) = 0;
% swe_vic_L15_2 = double(swe_vic_L15.*snsr_mask(:));

nt = size(swe_vic_L15, 2);
temp_outdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/swe_livneh';
mkdir(temp_outdir);
for t=1:nt
    swe_vic_map = fliplr(xyz2grid(info.lon, info.lat, swe_vic_L15(:,t)));
    geotiffwrite(fullfile(temp_outdir, ['swe_vic_' num2str(t) '.tif']), swe_vic_map, R);
end

swe_vic_map_avg_L15 = nanmean(swe_vic_L15, 2);
swe_vic_map_avg_L15 = fliplr(xyz2grid(info_L15.lon, info_L15.lat, swe_vic_map_avg_L15));
swe_vic_avg_ts_L15 = nanmean(swe_vic_L15,1)';

figure, imagesc(swe_vic_map_avg_L15)

%% Load SNSR data

fnames = dir('/Volumes/HD3/SNSR/cropped_resampled_tuo/*_cr.tif');

% Read data from WY1985 DOWY1 (10/1/1984) to WY2011 DOWY_last (9/30/2011)
swe_snsr = zeros(11,32,9861);
k = 1;
for wy_num = 1985:2011
    
    fnames1 = dir(['/Volumes/HD3/SNSR/cropped_resampled_tuo/*WY' num2str(wy_num) '*_cr.tif']);
    ndays_in_wy = length(fnames1);
    
    disp(['Processing WY ' num2str(wy_num)])
    for dowy_num = 1:ndays_in_wy
        fname = ['SN_SWE_WY' num2str(wy_num) '_dowy_' num2str(dowy_num) '_cr.tif'];
        tmp2 = geotiffread(fullfile('/Volumes/HD3/SNSR/cropped_resampled_tuo/', fname));
        tmp2(tmp2 < 0) = NaN;
        swe_snsr(:,:,k) = flipud(tmp2);
        k = k + 1;
    end
    
end
    
% nlat = length(lat);
% nlon = length(lon);
% swe_snsr = NaN(nlat, nlon, ndays);
% for k=1736:(ndays-91)
%     fname = ['SN_SWE_WY' num2str(WY(k)) '_dowy_' num2str(DOWY(k)) '_cr.tif'];
%     tmp2 = geotiffread(fullfile('/Volumes/HD3/SNSR/cropped_resampled_tuo/', fname));
%     tmp2(tmp2 < 0) = NaN;
%     swe_snsr(:,:,k) = flipud(tmp2);
%     disp(round(100*k/ndays, 2))   
% end

% average time series
% swe_ts_snsr = squeeze(nanmean(nanmean(swe_snsr, 1),2));

swe_ts_snsr = zeros(9861,1);
for k=1:9861
    tmp1 = swe_snsr(:,:,k);
    swe_ts_snsr(k) = nanmean(tmp1(:));
end

dlmwrite('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/overlapping_times/Final/swe_ts_snsr.txt', swe_ts_snsr, 'delimiter', ',')

% average map
swe_map_snsr = nanmean(swe_snsr,3);
figure, imagesc(swe_map_snsr)
geotiffwrite('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/overlapping_times/Final/swe_map_snsr.tif', flipud(swe_map_snsr), R)

% double check mean values
mean(swe_ts_snsr)
nanmean(swe_map_snsr(:))

% Re-calculate differently



% For some reason, the average map does not match the average time series
% Figure out why and recalculate. This problem does not occur in R.

% Minimal reproducible example:
sample_swe = zeros(2,2,2);
sample_swe(:,:,1) = [4,3;2,NaN];
sample_swe(:,:,2) = [0,7;0,NaN];

sample_swe_ts = squeeze(nanmean(nanmean(sample_swe, 1),2));
sample_swe_map = nanmean(sample_swe,3);

sample_swe_ts = zeros(2,1);
for k=1:2
    tmp1 = sample_swe(:,:,k);
    sample_swe_ts(k) = nanmean(tmp1(:));
end

%% Calculate SNSR time series

% Only take dates that overlap with the VIC simulation

[WY, DOWY] = CY2WY(vic_timevector);
ndays = length(vic_timevector);

nlat = length(lat);
nlon = length(lon);
swe_snsr = NaN(nlat, nlon, ndays);
for k=1736:(ndays-91)
    fname = ['SN_SWE_WY' num2str(WY(k)) '_dowy_' num2str(DOWY(k)) '_cr.tif'];
    tmp2 = geotiffread(fullfile('/Volumes/HD3/SNSR/cropped_resampled_tuo/', fname));
    tmp2(tmp2 < 0) = NaN;
    swe_snsr(:,:,k) = flipud(tmp2);
    disp(round(100*k/ndays, 2))   
end

save('/Volumes/HD3/SNSR/tuo_swe_snsr.mat', 'swe_snsr')

swe_map_snsr = nanmean(swe_snsr, 3);
swe_ts_snsr = squeeze(nanmean(nanmean(swe_snsr, 1)));

% figure, plotraster(lon, lat, swe_snsr(:,:,1952), 'SNSR SWE (mm)', '','')
% 
% swe_vic_map = fliplr(xyz2grid(info.lon, info.lat, swe_vic(:,1952)));
% figure
% plotraster(info.lon, info.lat, swe_vic_map, 'VIC SWE','Lon','Lat')

figure, plotraster(lon, lat, swe_map_snsr, 'Average SNSR SWE (mm)', '','')
figure, jsplot(vic_timevector, swe_ts_snsr, 'Average SNSR SWE', 'Time', 'SWE (mm)', 18)


%% Do the comparison

diffts = swe_ts_snsr - swe_vic_avg_ts;
diffmap = swe_map_snsr - swe_vic_map_avg;

%% Figure for ESSD paper

snsr_mask(snsr_mask==0) = NaN;
ts_nan_ind = isnan(swe_ts_snsr);

vic.time = vic_timevector(~ts_nan_ind);
vic.livneh.swe_ts = swe_vic_avg_ts_L15(~ts_nan_ind);
vic.vg.swe_ts = swe_vic_avg_ts(~ts_nan_ind);
snsr.swe_ts = swe_ts_snsr(~ts_nan_ind);
snsr.time = vic.time;

% Recalculate annual average SWE maps, account for the correct time range
vic.vg.swe = swe_vic(:,~ts_nan_ind);
vic.livneh.swe = swe_vic_L15(:,~ts_nan_ind);

vic.vg.avg_swe = nanmean(vic.vg.swe, 2);
vic.vg.avg_swe = fliplr(xyz2grid(info.lon, info.lat, vic.vg.avg_swe));

vic.livneh.avg_swe = nanmean(vic.livneh.swe, 2);
vic.livneh.avg_swe = fliplr(xyz2grid(info_L15.lon, info_L15.lat, vic.livneh.avg_swe));

SWE1 = flipud(vic.livneh.avg_swe).*snsr_mask; % L2013
SWE2 = flipud(vic.vg.avg_swe).*snsr_mask; % VICGlobal
SWE3 = swe_map_snsr.*snsr_mask; % SNSR

%% Make figure for paper (final)

indir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/overlapping_times/Final';

swe_ts_livneh = load(fullfile(indir, 'swe_ts_livneh.txt'));
swe_ts_vicglobal = load(fullfile(indir, 'swe_ts_vicglobal.txt'));
% swe_ts_snsr

swe_map_livneh = flipud(geotiffread(fullfile(indir, 'swe_livneh.tif')));
swe_map_livneh(swe_map_livneh<0) = NaN;
swe_map_vicglobal = flipud(geotiffread(fullfile(indir, 'swe_vicglobal.tif')));
swe_map_vicglobal(swe_map_vicglobal<0) = NaN;
% swe_map_snsr

timevector = datetime(1984, 10, 1):datetime(2011,9,30);
%%
figure

subplot(2,3,1)
plotraster(lon, lat, swe_map_livneh, 'L2013 SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,2)
plotraster(lon, lat, swe_map_vicglobal, 'VICGlobal SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,3)
plotraster(lon, lat, swe_map_snsr, 'SNSR SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,[4 6])
linewidth = 1.5;
hold on
plot(timevector, swe_ts_livneh, 'linewidth', linewidth, 'color', 'blue')
plot(timevector, swe_ts_vicglobal, 'linewidth', linewidth, 'color', 'red')
plot(timevector, swe_ts_snsr, 'linewidth', linewidth, 'color', 'black')
xlabel('Time')
ylabel('SWE (mm)')
% title('SWE comparison')
set(gca, 'fontsize', 20)
grid on
legend('L2013', 'VICGlobal', 'SNSR', 'Location', 'North', 'Orientation', 'Horizontal')


%% Calculate RMSE

addpath('/Volumes/HD3/SWOTDA/Calibration/CalVIC')
RMSE.vg = myRMSE(snsr.swe_ts, vic.vg.swe_ts);
RMSE.livneh = myRMSE(snsr.swe_ts, vic.livneh.swe_ts);

%% Save results of calculations

mkdir('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/')
save('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/swe_data.mat')
% clear
% load('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/SWE_Validation/swe_data.mat')

%%
figure, 
subplot(2,1,2)
jsplot(vic_timevector, diffts, 'SNSR - VICGlobal', 'Time', 'SWE (mm)', 18)
ylim([0,300])
grid on

figure
subplot(3,1,1)
plotraster(lon, lat, swe_map_snsr, 'Average SNSR SWE (mm)', 'Lon','Lat')
% caxis([0, 500])
subplot(3,1,2)
% swe_vic_map_avg_copy = swe_vic_map_avg;
% swe_vic_map_avg_copy(swe_vic_map_avg<1) = NaN;
plotraster(lon, lat, flipud(swe_vic_map_avg), 'Average VIC SWE (VICGlobal, mm)', 'Lon','Lat')
% caxis([0, 500])
% subplot(4,1,3)
% plotraster(lon, lat, flipud(swe_vic_map_avg_L15), 'Average VIC SWE (Livneh, mm)', 'Lon','Lat')
subplot(3,1,3)
plotraster(lon, lat, diffmap, 'Difference SNSR - VICGlobal SWE (mm)', 'Lon','Lat')
% caxis([-250,250])
% colormap(bluewhitered)
% caxis([0, 500])

%% Figure for Dongyue

figure
linewidth = 2;
plot(vic.time, vic.livneh.swe_ts, 'linewidth', linewidth, 'color', 'blue')
hold on
plot(vic.time, vic.vg.swe_ts, 'linewidth', linewidth, 'color', 'red')
xlabel('Time')
ylabel('SWE (mm)')
% title('SWE comparison')
set(gca, 'fontsize', 18)
grid on
legend('L2013', 'VICGlobal', 'Location', 'North', 'Orientation', 'Horizontal')

%% Repeat figures with the updated data

saveloc = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation';
[swe_livneh, R2] = geotiffread(fullfile(saveloc, 'swe_livneh.tif'));
swe_vicglobal = geotiffread(fullfile(saveloc, 'swe_vicglobal.tif'));

swe_livneh(swe_livneh<0) = NaN;
swe_vicglobal(swe_vicglobal<0) = NaN;

swe_livneh = flipud(swe_livneh);
swe_vicglobal = flipud(swe_vicglobal);

swe_livneh_ts = load(fullfile(saveloc, 'swe_ts_livneh.txt'));
swe_vicglobal_ts = load(fullfile(saveloc, 'swe_ts_vicglobal.txt'));

%% Make the figure

figure

subplot(2,3,1)
plotraster(lon, lat, swe_livneh, 'L2013 SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,2)
plotraster(lon, lat, swe_vicglobal, 'VICGlobal SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,3)
plotraster(lon, lat, SWE3, 'SNSR SWE (mm)', 'Lon','Lat')
caxis([0, 530])

subplot(2,3,[4 6])
linewidth = 1;
hold on
plot(vic_timevector, swe_livneh_ts, 'linewidth', linewidth, 'color', 'blue')
plot(vic_timevector, swe_vicglobal_ts, 'linewidth', linewidth, 'color', 'red')
plot(vic.time, snsr.swe_ts, 'linewidth', linewidth, 'color', 'black')
xlabel('Time')
ylabel('SWE (mm)')
% title('SWE comparison')
set(gca, 'fontsize', 18)
grid on
legend('L2013', 'VICGlobal', 'SNSR', 'Location', 'North', 'Orientation', 'Horizontal')
