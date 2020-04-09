% Makes plots of soil parameters from the VIC4 soil parameter file
%
% tuolumne inputs
% tifdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_VICGlobal/soil_plots/tifs';
% outdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_VICGlobal/soil_plots/pngs';
%
% one-degree tile inputs
% tifdir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/soil_plots/tifs';
% outdir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/soil_plots/pngs';

function plot_soils(tifdir, outdir)

tifnames = dir([tifdir '/*.tif']);

[var, R] = geotiffread(fullfile(tifdir, tifnames(1).name));
Rmat = georefobj2mat(R);
[lon, lat] = pixcenters(Rmat, size(var));

% [elev, R] = geotiffread(fullfile(tifdir, tifnames(17).name));
% figure, plotraster(lon, lat, elev, 'elev','', '')
% elev(1,1)

for j=1:length(tifnames)
    var = flipud(geotiffread(fullfile(tifdir, tifnames(j).name)));
    varname = strsplit(tifnames(j).name, '.');
    varname = varname{1};
    f1 = figure;
    plotraster(lon, lat, var, varname, 'Lon', 'Lat');
    saveas(f1, fullfile(outdir, [varname '.png']))
    close(f1)
end

return