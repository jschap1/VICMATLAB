% Script for snapping soil parameter file to the coordinates of forcing
% data files, which may be on a slightly different grid


soils_clipped = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_L2015/soils_tuolumne.txt'); % load soil parameter file

[lon, lat] = get_coordinates_from_VIC_file('/Volumes/HD4/SWOTDA/Data/Tuolumne/forc_ascii', 'Forcing');

latlon = [lat, lon];
precision = 5;
outformat = 'livneh';
outname = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_L2015/soils_tuolumne_snapped.txt';
[soils_clipped, newcoords, oldcoords] = snap_soils(soils_clipped, latlon, precision, outname, outformat);