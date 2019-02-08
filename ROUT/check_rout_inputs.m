% June 25, 2018, JRS
% Double check routing model inputs

[fract, ~] = arcgridread('/Volumes/HD3/SWOTDA/Data/UMRB/ROUT/umrb.fract');
[fdir, ~] = arcgridread('/Volumes/HD3/SWOTDA/Data/UMRB/ROUT/umrb.flowdir');

stationsfile = '/Volumes/HD3/SWOTDA/Data/UMRB/ROUT/umrb.stations';
f1 = fopen(stationsfile);
stations = textscan(f1, '%d %s %f %f %f');
fclose(f1);

% There should not be negative values in the stations file...

% Plot points where stations are
% Use VIC routing model convention

stations{3} % col, from left
stations{4} % row, from bottom

[nrow, ncol] = size(fract);

% Plot flow direction
figure, imagesc(fdir), colorbar, title('flowdir')

colfromleft = stations{3};
rowfrombottom = stations{4};
[xx,yy] = pix2latlon(R, colfromleft, rowfrombottom)


%%

fdir = geotiffread('/Volumes/HD3/SWOTDA/Data/UMRB/Flowdir/umrb_flowdir.tif');
fdir = arcgridread('/Volumes/HD3/SWOTDA/Data/UMRB/ROUT/umrb.flowdir');
fract = arcgridread('/Volumes/HD3/SWOTDA/Data/UMRB/ROUT/umrb.fract');
fdir_grass = convertflowdir(fdir, 'grass');

imagesc(fdir);
imagesc(fdir_grass);