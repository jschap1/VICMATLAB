%% Remove negative values from the classic mode snowband file

snowband_classic = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/snowbands_MERIT.txt');
snowband_classic_old = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/Old/snowbands_MERIT.txt');
minvals = zeros(16,1);
for k=1:16
    minvals(k) = min(snowband_classic_old(:,k));
end
ind = find(snowband_classic_old(:,7:11)<0);
snowband_classic_old(snowband_classic_old<0) = 0;

% Write out the snowband file
outfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/snowbands_MERIT_min0.txt';
dlmwrite(outfile, snowband_classic_old, 'Delimiter', '\t','precision', 8);

%% Remove negative values from the image mode parameter file

ncname = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Image/VICGlobal_params.nc';
snowband_image = ncread(ncname, 'elevation');
lat = ncread(ncname, 'lat');
lon = ncread(ncname, 'lon');

ind = find(snowband_image<0);
snowband_image(snowband_image<0) = 0;

% Write out the parameter file with the corrected snowband data

ncwrite(ncname, 'elevation', snowband_image);