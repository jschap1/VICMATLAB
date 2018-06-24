% June 25, 2018, JRS
% Double check routing model inputs

[fract, ~] = arcgridread('/Users/jschap/Desktop/UMRB/R_In/umrb.fract');
[fdir, ~] = arcgridread('/Users/jschap/Desktop/UMRB/R_In/umrb.flowdir');

stationsfile = '/Users/jschap/Desktop/UMRB/R_In/umrb.stations_noNONE.txt';
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
