% Prepare MERRA-2 Forcing Data
%
% Prepares forcing data for the VIC model using NetCDF MERRA-2 data
% Downscales MERRA-2 data and crops it to the study domain
%
% INPUTS
% MERRA-2 atmospheric reanalysis data
% Specifically, flx, rad, and slv files
% Fluxes: M2T1NXFLX: MERRA-2 tavg1_2d_flx_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics V5.12.4
% Radiation: M2T1NXRAD: MERRA-2 tavg1_2d_rad_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics V5.12.4
% Single Level Diagnostics: M2T1NXSLV: MERRA-2 tavg1_2d_slv_Nx: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4
%
% Updated 3/16/2019 JRS, GK
% d3 - added cropping functionality so it only writes out data for the basin extent
% d5 - double checked units, added comments, modified conversion formulas
% Updated 4/15/2019 JRS to perform cropping before interpolating to
% speed up runtime.
% Updated 4/24/2019 v4 JRS to dramatically increase speed by changing the way
% it writes out the forcing input files
%
% Updated 4/24/2019 v5
% To reduce RAM-pressure, it runs over one variable at a time
%
% Updated 4/25/2019 v6 
% Uses tall arrays to avoid having to load in all the data at once
%
% Updated 4/29/2019 v67
% Using plain old datastores, instead of tall arrays. This seems the best
% solution.

%% Set working directory

% user = 'Gurjot Kohli';
user = 'jschap';

if strcmp(user, 'jschap')
%     cd('/Volumes/HD_ExFAT/MERRA2')
    cd('/Users/jschap/Desktop/MERRA2/Raw/Sample2')
    addpath('/Users/jschap/Desktop/MERRA2/Codes')
    addpath('/Users/jschap/Documents/Codes/VICMATLAB')
elseif strcmp(user, 'Gurjot Kohli')
    cd('C:\Users\Gurjot Kohli\Box\Shared_SWOTDA')
elseif strcmp(user, 'elqui')
    cd('/Volumes/elqui_hd5/nobackup/DATA/MERRA2');
    addpath('/Volumes/LIT')
end

%% Inputs

delete(gcp('nocreate')) % remove any existing parallel pools
% p = parpool(); % start a parallel pool

lat_range = [24 37.125];
lon_range = [66 82.875];

% lat_range = [24 25];
% lon_range = [66 67];

% lat_range = [24.02500 37.08333];
% lon_range = [66.15833 82.45000];

% lat_range = [-90 90];
% lon_range = [-180 180];

target_res = 1/16; % can change this if not enough RAM

% MERRA-2 static file with lat/lon info
static_file = '/Users/jschap/Desktop/MERRA2/Processed/static_file_MERRA2.in'; 

% Directory to save the downscaled forcing files
forcdir = '/Users/jschap/Desktop/MERRA2/Forc_1980-2019'; 

% MERRA-2 filenames
prec_names = dir('MERRA2*flx*.hdf');
rad_names = dir('MERRA2*rad*.hdf');
air_names = dir('MERRA2*slv*.hdf');

ndays = length(prec_names);
% ndays = 10;

% soils = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/soils_clipped.txt');
% slat = soils(:,3);
% slon = soils(:,4);

xres = 5/8;
yres = 1/2;

%% load static data
A = load(static_file);
lonlat = [A(:,3), A(:,2)];
lon = sort(unique(lonlat(:,1)));
lat = sort(unique(lonlat(:,2)));

%% Make cropping rectangle

% must use grid cell centers for makerefmat
R1 = makerefmat(min(lon), min(lat), xres, yres);

%     this formulation guarantees whole number indices
%     [ymin, xmin] = latlon2pix(R1, min(lat)+yres*10, min(lon)+xres*5) % lat, lon

% find appropriate minimum values for the cropping rectangle
minval_opt_1 = 10;
minval_opt_2 = 10;
for p=1:400
    tmp1 = abs(min(lon)+ p*xres - lon_range(1));
    tmp2 = abs(min(lat)+ p*yres - lat_range(1));
    if tmp1 < minval_opt_1
        minval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < minval_opt_2
        minval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
minlon = min(lon)+ p_opt_1*xres;
minlat = min(lat)+ p_opt_2*yres;

% make a 1-pixel border just to be safe
minlon = minlon - xres;
minlat = minlat - yres;

% find appropriate maximum values for the cropping rectangle
maxval_opt_1 = 10;
maxval_opt_2 = 10;
for p=1:400
    tmp1 = abs(minlon + p*xres - lon_range(2));
    tmp2 = abs(minlat + p*yres - lat_range(2));
    if tmp1 < maxval_opt_1
        maxval_opt_1 = tmp1;
        p_opt_1 = p;
    end
    if tmp2 < maxval_opt_2
        maxval_opt_2 = tmp2;
        p_opt_2 = p;
    end
end
maxlon = minlon + p_opt_1*xres;
maxlat = minlat + p_opt_2*yres;  

% make a 1-pixel border just to be safe
maxlon = maxlon + xres;
maxlat = maxlat + yres;

[ymin, xmin] = latlon2pix(R1, minlat, minlon);
[ymax, xmax] = latlon2pix(R1, maxlat, maxlon);
width = xmax - xmin;
height = ymax - ymin;

%     rect = [floor(xmin), floor(ymin), ceil(width), ceil(height)];
rect = [xmin, ymin, width, height];

% I have to do this to get the code to work...
height = height+1;
width = width+1;

out_lat = minlat:yres:maxlat;
out_lon = minlon:xres:maxlon;
    
%% Downscale the meteorological forcings

prec_fine = cell(ndays, 1);
ps_fine = cell(ndays, 1);
swdn_fine = cell(ndays, 1);
lwdn_fine = cell(ndays, 1);
vp_fine = cell(ndays, 1);
wind_fine = cell(ndays, 1);

%% get datetime information for the MERRA-2 data
datetimearray = zeros(ndays, 4);
for d=1:ndays
    air_name = air_names(d).name;
    mn_split = strsplit(air_name, '.');
    mn_date = mn_split{3};
    yr = str2double(mn_date(1:4));
    mon = str2double(mn_date(5:6));
    day = str2double(mn_date(7:8));
    datetimearray(d,:) = [yr, mon, day, 0];
end
forc_dates = datetime(datetimearray(:,1), datetimearray(:,2), datetimearray(:,3));
forc_date_string = datestr(forc_dates); % useful for making filenames

%% Crop, downscale the forcing data (temperature)

% temporary location to store the downscaled, cropped forcing data as text files
% intermediate_dir = '/Volumes/HD_ExFAT/downscaled_cropped_forcings';
intermediate_dir = '/Users/jschap/Desktop/MERRA2/Downscaled';
merra_dir = '/Users/jschap/Desktop/MERRA2/Raw/Sample2';

return_cell = 0;
for d=1:ndays
    
    temperature = hdfread(fullfile(merra_dir, air_names(d).name), 'EOSGRID', 'Fields', 'T2M');

    temp_crop = zeros(24, height, width);
    for h=1:24
        temp_crop(h,:,:) = imcrop(squeeze(temperature(h,:,:)), rect);
    end

    if d==1
        [temp_fine, target_lon, target_lat] = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
        [nr, nc, ~] = size(temp_fine);
    else
        temp_fine = interpolate_merra2(temp_crop, target_res, out_lon, out_lat, return_cell);
    end
    
    % write out temp_fine arrays for each day, then read them back in as a
    % datastore in order to manipulate them and write VIC input files
    
    temp_fine_mat = reshape(temp_fine, nr*nc, 24); % converting to 2D
    savename = ['temperature_', forc_date_string(d,:), '.txt'];
%     dlmwrite(fullfile(intermediate_dir, savename), temp_fine_mat') % nt by ncells
    dlmwrite(fullfile(intermediate_dir, savename), temp_fine_mat) % ncells by nt
    
%     temp_fine2 = reshape(temp_fine_mat, 240, 304, 24); % converting back to 3D
    
end

%% Load the cropped, downscaled forcing data and prepare VIC input files

ds = tabularTextDatastore(fullfile(intermediate_dir, 'temperature*'), 'FileExtension', '.txt');

mapreducer(0); % do not use parallel pool
% took 29 minutes with 6 processors for ndays = 2
% took 22 minutes with one processor (lots of overhead time for // pool)

% make list of 72960 latlon values corresponding to each grid cell, and in
% the order they appear in each of the saved temperature.txt files
xyz = raster2xyz(target_lon', target_lat', ones(nr, nc));
lonlat = xyz(:,1:2);
% hopefully the order is correct; plotting will show whether it worked

ncells = nr*nc; % number of cells in the cropped, downscaled forcing data
precision = '%3.5f';

% forcings_out = cell(ncells, 1);
% mapreducer(gcp);

tic; % this step takes a while, on the order of an hour
ta = tall(ds); % the more columns there are (ncells), the longer this step takes
toc % took 400 seconds for the full IRB domain and ndays=5
% p.IdleTimeout = 10; % go to // preferences to set

% If necessary for computational considerations, change the shape of the tall array to npix by nt and
% use fewer timesteps. Disaggregate to 3-hourly and use a shorter modeling
% period, instead of the full record.

% Make forcing file for one grid cell
tic; 
forcing_out = zeros(ndays*24,7);
reset(ds);
ds.ReadSize = 'file';
for d=1:ndays
    aa = read(ds);
    d1 = 24*(d-1)+1;
    d2 = 24*d;    
    forcing_out(d1:d2,1) = table2array(aa(1,:))';
end
toc
% 1.56 seconds
% BOOM. Done.
% Over more days, this will take a bit longer, but should still be
% reasonable.

% Do for all grid cells


% calculate indices ahead of time
d1 = zeros(ndays,1);
d2 = zeros(ndays,1);
for d=1:ndays
    d1(d) = 24*(d-1)+1;
    d2(d) = 24*d;   
end
        
% This will take a couple of days for the full domain
% A parfor loop might be helpful for cutting execution time 
tic; 
% parpool()
ds.ReadSize = 'file';
parfor x=1:ncells
    savename1 = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) '.txt']);
    forcing_out = zeros(ndays*24,7);
    reset(ds);
    for d=1:ndays
        aa = read(ds);
        forcing_out(d1(d):d2(d),1) = table2array(aa(x,:))'; % temperature
    end
    dlmwrite(savename1, forcing_out)
end
toc;
% 148.5 s for 100 cells
% 30.5 hours for 73000 cells

forcing_out = zeros(ndays*24,7);
reset(ds);
ds.ReadSize = 'file';
for d=1:ndays
    aa = read(ds);
    d1 = 24*(d-1)+1;
    d2 = 24*d;    
    forcing_out(d1:d2,1) = table2array(aa(1,:))';
end
toc



tic
for x=1:ncells

    savename = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) '.txt']);
%     forcing_temp = ta(:,x);

%     aa = gather(ta(:,x)); % takes about 4.3 seconds with 120 hours

    tmp = gather(ta(:,x));
    forcings_out = table2array(tmp); % about 6s per call with 2304 cells by 120 hours
    dlmwrite(savename, forcings_out)
    
    if mod(x,100)==0
        disp(x)
    end
    
%     forcings_out = zeros(24*ndays, 7);
%     forcings_out(1,:) = gather(ta(:,x)); % took about 12 minutes on one processor for ncells=70e3
   
%     forcings_out(1,:) = ta(:,x); % try using gather outside the loop (it will take too long otherwise)
%     mapreducer(0) w/// toolbox, this seems to take a while (30 min.). Profile it.
%     forcings_out(1,:) = ta(:,x);

    % change this to have all 7 different forcings
    
%     write(fullfile(intermediate_dir, 'sdf'), forcings_out)
    
end
toc
% aa = gather(forcings_out);
% if this takes 3s per cell, then it takes about 2 hours for 2304 cells
% calling gather within a loop really slows things down
% Using more workers in a parallel pool speeds up how long it takes to call
% gather. Takes under 1s per grid cell with 6 workers (40 min for stack)
%
% If this continues to take too long, try running on Hoffman or another
% parallel computing resource, such as Matlab Parallel Server or AWS.
% Hoffman2 has Matlab R2015b. Tall arrays were intro in 2016b.
%
% As long as this takes less than __ seconds per grid cell, it is
% acceptable

    % get temperature for a grid cell all the way across the array  
    for d=1:ndays
        temp_day = zeros(24, 1);
        for h=1:24
            temp_day(h) = ta();
            temp_day(h) = temp_fine{d}{h}(row1, col1); % avoid reading the datastore ndays times
        end
        temp_day = temp_day - 273.15; % Kelvin to Celsius
        d1 = 24*(d-1)+1;
        d2 = 24*d;
        forcings_out(d1:d2,:) = [temp_day, prec_day, ps_day, swdn_day, lwdn_day, vp_day, wind_day];        
    end
    

for x=1:nx % loop over grid cells
    for y=1:ny 
    
        for k=1:ndays   % loop over days      
            [row1, col1] = latlon2pix(R, target_lat(y), target_lon(x));
            temp_day = zeros(24, 1);
            for h=1:24
                temp_day(h) = temp_fine{k}{h}(row1, col1);
            end
            temp_day = temp_day - 273.15; % Kelvin to Celsius
            
            d1 = 24*(k-1)+1;
            d2 = 24*k;
            forcings_out(d1:d2,:) = [temp_day, prec_day, ps_day, swdn_day, lwdn_day, vp_day, wind_day];
        end
        % write forcing data for this gridcell to file
        fID = fopen(savename, 'w');
        formatstring = '%0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f\n';
        fprintf(fID, formatstring, forcings_out');
        fclose(fID);     
    end
end



% fcn = @(fname) hdfread(fname, 'EOSGRID', 'Fields', 'T2M');

% Good resource: https://slideplayer.com/slide/12846850/

% fds = fileDatastore('/Volumes/HD_ExFAT/Sample', 'ReadFcn', @h5readAll);
% Turn off the parallel computing option for now
% See https://www.mathworks.com/help/parallel-computing/run-tall-arrays-on-a-parallel-pool.html
mapreducer(0);

% create a tall array of the temperature data
merra_dir = '/Volumes/HD_ExFAT/Sample';
merra_dir = '/Volumes/HD_ExFAT/MERRA2';
convert_hdf_to_mat(merra_dir) % I don't think this is necessary. Can just read HDF files directly
fds = fileDatastore(merra_dir, 'ReadFcn', @load, 'FileExtensions', '.mat');

temp_fine = zeros(240, 304, 24); % need to get 240, 304 programmatically
% cannot actually make temp_fine the full size (24*ndays); it would be too large
% combine multiples?

reset(fds)
temperature_tall = tall(fds);


a = readall(fds);
a{1}


ta2 = cell2underlying(ta);

ta = tall(fds);
ta(1) % data for day 1
ta(2) % data for day 2

m = gather(mean(ta));

% crop each one
nhours = 24;
temp_crop = zeros(nhours, height, width);
for h=1:nhours
    temp_crop(h,:,:) = imcrop(squeeze(ta(1)(h,:,:)), rect);
end
% note: tall cell arrays are less useful than tall tabular arrays


% crop the data
fds.ReadSize = 'file';
raw_merra = read(fds);
raw_merra = raw_merra.dat1;



[temp_fine{k}, target_lon, target_lat] = interpolate_merra2(temp_crop, target_res, out_lon, out_lat);


temp_fine = cell(ndays, 1); % quickly grows to 10s of GB... Not OK. 

for k=1:ndays
        
    air_name = air_names(k).name;
    temp = hdfread(air_name, 'EOSGRID', 'Fields', 'T2M');
    
    nhours = size(temp,1);
    temp_crop = zeros(nhours, height, width);
    for h=1:nhours
        temp_crop(h,:,:) = imcrop(squeeze(temp(h,:,:)), rect);
    end
    
    [temp_fine{k}, target_lon, target_lat] = interpolate_merra2(temp_crop, target_res, out_lon, out_lat);

    nx = length(target_lon);
    ny = length(target_lat);
    R = makerefmat(target_lon(1), target_lat(1), target_res, target_res);
    
    % Write dates to file
    % This is a check to make sure the VIC forcing files are in the
    % right order and aren't missing any data
    fID2 = fopen(fullfile(forcdir, 'forcing_dates.txt'), 'a');
    dates_out = [str2double(yr), str2double(mon), str2double(day)];
    formatstring = '%d %d %d\n';
    fprintf(fID2, formatstring, dates_out');
    fclose(fID2);
    
    disp(['Interpolated data for day ' num2str(k) ' of ' num2str(ndays)]) % progress tracker
    
    % Compile files for each pixel across times
    % Write data to files for each pixel
       
end

save('temperature_fine_interp.mat', temp_fine, '-v7.3')
save('temperature_fine_interp.mat', temp_fine)
  

%%
for k=1:ndays % do in parallel for individual MERRA-2 files?
    
    % get information about the MERRA-2 file
    prec_name = prec_names(k).name;
    mn_split = strsplit(prec_name, '.');
    mn_date = mn_split{3};
    yr = mn_date(1:4);
    mon = mn_date(5:6);
    day = mn_date(7:8);
        
    air_name = air_names(k).name;
    rad_name = rad_names(k).name;
        
    % load meteorological forcing variables   
    prec = hdfread(prec_name, 'EOSGRID', 'Fields', 'PRECTOT');
    
    temp = hdfread(air_name, 'EOSGRID', 'Fields', 'T2M');
    
    ps = hdfread(air_name, 'EOSGRID', 'Fields', 'PS'); % pressure
    
    qv = hdfread(air_name, 'EOSGRID', 'Fields', 'QV2M'); % specific humidity
%     qv2vp = @(qv,ps) (qv.*ps)./((0.378.*qv)+0.622);
    qv2vp = @(qv,ps) (qv.*ps)./(qv+0.622); % vp units are same as ps units
    vp = qv2vp(qv, ps); % convert to vapor pressure
    
    windx = hdfread(air_name, 'EOSGRID', 'Fields', 'U2M');
    windy = hdfread(air_name, 'EOSGRID', 'Fields', 'V2M');
    get_resultant = @(x,y) (x.^2 + y.^2).^(1/2);
    wind = get_resultant(windx, windy); % get resultant wind vector
    
    swdn = hdfread(rad_name, 'EOSGRID', 'Fields', 'SWGDN');
    
    % Long wave calculation (Idso, 1981; overly simplistic)

    eo = vp/100; % surface level vapor pressure (mb), converting from Pascals to millibars 
%     eo = ps.*qv./0.622; 
    sigma_sb = 5.670e-8; % Stephen-Boltzmann constant
    emissivity =  0.179*(eo./100).^(1/7).*(exp(350./temp)); % temperature (K)
    lwdn = emissivity.*sigma_sb.*temp.^4; % (W/m^2)
    % The major issue here is that the formula assumes a cloudless atmosphere
    
    %% Crop to modeling extent
         
    % crop each variable
    nhours = size(prec,1);
    
    prec_crop = zeros(nhours, height, width);
    temp_crop = zeros(nhours, height, width);
    ps_crop = zeros(nhours, height, width);
    vp_crop = zeros(nhours, height, width);
    wind_crop = zeros(nhours, height, width);
    swdn_crop = zeros(nhours, height, width);
    lwdn_crop = zeros(nhours, height, width);
    for h=1:nhours
        prec_crop(h,:,:) = imcrop(squeeze(prec(h,:,:)), rect);
        temp_crop(h,:,:) = imcrop(squeeze(temp(h,:,:)), rect);
        ps_crop(h,:,:) = imcrop(squeeze(ps(h,:,:)), rect);
        vp_crop(h,:,:) = imcrop(squeeze(vp(h,:,:)), rect);
        wind_crop(h,:,:) = imcrop(squeeze(wind(h,:,:)), rect);
        swdn_crop(h,:,:) = imcrop(squeeze(swdn(h,:,:)), rect);
        lwdn_crop(h,:,:) = imcrop(squeeze(lwdn(h,:,:)), rect);
    end
    %     clear prec_fine temp_fine ps_fine vp_fine wind_fine swdn_fine lwdn_fine % to conserve RAM
    
%     out_lat = lat_range(1):1/2:lat_range(2);
%     out_lon = lon_range(1):5/8:lon_range(2);
    
    out_lat = minlat:1/2:maxlat;
    out_lon = minlon:5/8:maxlon;
        
%     figure, imagesc(lon_range, lat_range, swdn_crop{12})
%     colorbar, xlabel('Lon'); ylabel('Lat');
%     set(gca, 'ydir', 'normal')

% if doing the whole domain (global), then no need to crop
% prec_crop = prec;
% temp_crop = temp;
% ps_crop = ps;
% vp_crop = vp;
% wind_crop = wind;
% swdn_crop = swdn;
% lwdn_crop = lwdn;
% out_lat = lat;
% out_lon = lon;

    %% Interpolate to target resolution
    
    % this can use a lot of RAM depending on the target_resolution
    [prec_fine{k}, target_lon, target_lat] = interpolate_merra2(prec_crop, target_res, out_lon, out_lat);
    temp_fine{k} = interpolate_merra2(temp_crop, target_res, out_lon, out_lat);
    ps_fine{k} = interpolate_merra2(ps_crop, target_res, out_lon, out_lat);
    vp_fine{k} = interpolate_merra2(vp_crop, target_res, out_lon, out_lat);
    wind_fine{k} = interpolate_merra2(wind_crop, target_res, out_lon, out_lat);
    swdn_fine{k} = interpolate_merra2(swdn_crop, target_res, out_lon, out_lat);
    lwdn_fine{k} = interpolate_merra2(lwdn_crop, target_res, out_lon, out_lat);
    
    % Write out a GeoTIFF (optional)   
    nx = length(target_lon);
    ny = length(target_lat);
    R = makerefmat(target_lon(1), target_lat(1), target_res, target_res);
%     geotiffwrite('/Users/jschap/prec_fine_crop.tif', prec_fine{1}, R)

    % Write dates to file
    % This is a check to make sure the VIC forcing files are in the
    % right order and aren't missing any data

    fID2 = fopen(fullfile(forcdir, 'forcing_dates.txt'), 'a');
    dates_out = [str2double(yr), str2double(mon), str2double(day)];
    formatstring = '%d %d %d\n';
    fprintf(fID2, formatstring, dates_out');
    fclose(fID2);
    
    disp(['Interpolated data for day ' num2str(k) ' of ' num2str(ndays)]) % progress tracker
    
end
  
    
%% Write VIC input files

for x=1:nx % loop over grid cells
    for y=1:ny 
        
        % make savename for this gridcell
        precision = '%3.5f';
        savename = fullfile(forcdir, ['Forcings_' num2str(target_lat(y), precision) '_' num2str(target_lon(x), precision) '.txt']);        
        
        % initialize output forcing data
        % forcings_out = zeros(24*ndays, 7);

        for k=1:ndays   % loop over days      

            [row1, col1] = latlon2pix(R, target_lat(y), target_lon(x));

            prec_day = zeros(24, 1);
            temp_day = zeros(24, 1);
            ps_day = zeros(24, 1);
            vp_day = zeros(24, 1);
            wind_day = zeros(24, 1);
            swdn_day = zeros(24, 1);
            lwdn_day = zeros(24, 1);
            for h=1:24
                prec_day(h) = prec_fine{k}{h}(row1, col1);
                temp_day(h) = temp_fine{k}{h}(row1, col1);
                ps_day(h) = ps_fine{k}{h}(row1, col1);
                vp_day(h) = vp_fine{k}{h}(row1, col1);
                wind_day(h) = wind_fine{k}{h}(row1, col1);
                swdn_day(h) = swdn_fine{k}{h}(row1, col1);
                lwdn_day(h) = lwdn_fine{k}{h}(row1, col1);
            end

            % MERRA-2 original units
            % precipitation (kg m-2 s-1)
            % temperature (K)
            % pressure (Pa)
            % downwelling shortwave radiation (W/m^2)
            % wind (m/s)
            % specific humidity (kg/kg)

            % Units of calculated variables
            % vapor pressure (???)

            % downwelling longwave radiation (W/m^2)

            % Unit conversions for VIC
            temp_day = temp_day - 273.15; % Kelvin to Celsius
            ps_day = ps_day/1000; % Pascal to kPa
            vp_day = vp_day/1000; % Pascal to kPa
            prec_day = 3600*prec_day; % kg/m2/s to mm/hr

            % units should be:
            % precipitation (mm/timestep)
            % temperature (deg. C)
            % pressure (kPa)
            % downwelling shortwave radiation (W/m^2)
            % downwelling longwave radiation (W/m^2)
            % vapor pressure (kPa)
            % wind (m/s)
            
            d1 = 24*(k-1)+1;
            d2 = 24*k;
            forcings_out(d1:d2,:) = [temp_day, prec_day, ps_day, swdn_day, lwdn_day, vp_day, wind_day];

        end
        
        % write forcing data for this gridcell to file
        fID = fopen(savename, 'w');
        formatstring = '%0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f\n';
        fprintf(fID, formatstring, forcings_out');
        fclose(fID);
            
    end
end
   

%% Plot forcing time series for a particular grid cell

% The forcing files have either 1920 or 1944 entries (80 or 81 days)
% They start on Oct. 2, 2017 and go until Dec. 20 or 21, 2017

dat = load('/Users/jschap/Desktop/MERRA2/Forc/Forcings_25.46875_70.46875.txt');
dat2 = load('/Users/jschap/Desktop/MERRA2/Forc/Forcings_37.65625_83.40625.txt');

temperature = dat(:,1);
precipitation = dat(:,2);
pressure = dat(:,3);
shortwave = dat(:,4);
longwave = dat(:,5);
vapor = dat(:,6);
wind = dat(:,7);

figure, subplot(4,2,1)
plot(temperature), title('Temperature (degrees C)'), xlabel('Hours')

subplot(4,2,2)
plot(precipitation), title('Precipitation (mm/hr)'), xlabel('Hours')

nprecip = length(precipitation(precipitation>0.001));

subplot(4,2,3)
plot(pressure), title('Pressure (kPa)'), xlabel('Hours')

subplot(4,2,4)
plot(shortwave), title('Downwelling shortwave (W/m^2)'), xlabel('Hours')

subplot(4,2,5)
plot(longwave), title('Downwelling longwave (W/m^2)'), xlabel('Hours')

subplot(4,2,6)
plot(vapor), title('Vapor pressure (kPa)'), xlabel('Hours')

subplot(4,2,7)
plot(wind), title('Wind speed (m/s)'), xlabel('Hours')


%%


if (min(abs(slat - flat)) < 0.5*target_res) && (min(abs(slon - flon)) < 0.5*target_res)
    1;
end


% Loop through each grid cell and match to closest grid cell in soil
% parameter file


