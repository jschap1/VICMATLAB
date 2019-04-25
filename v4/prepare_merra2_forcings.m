% Prepare MERRA-2 Forcing Data
%
% Prepares forcing data for the VIC model using NetCDF MERRA-2 data
% Downscales MERRA-2 data and crops it to the study domain
%
% Updated 3/16/2019 JRS, GK
% d3 - added cropping functionality so it only writes out data for the basin extent
% d5 - double checked units, added comments, modified conversion formulas
% Updated 4/15/2019 JRS

%% Set working directory

% user = 'Gurjot Kohli';
user = 'jschap';

if strcmp(user, 'jschap')
    cd('/Volumes/HD_ExFAT/MERRA2')
    addpath('/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Codes')
elseif strcmp(user, 'Gurjot Kohli')
    cd('C:\Users\Gurjot Kohli\Box\Shared_SWOTDA')
elseif strcmp(user, 'elqui')
    cd('/Volumes/elqui_hd5/nobackup/DATA/MERRA2');
    addpath('/Volumes/LIT')
end

%% Downscale the meteorological forcings

% input extent to output data
% lat_range = [26.02500 28.08333]; 
% lon_range = [64.15833 66.45000];

lat_range = [24.02500 37.08333]; 
lon_range = [66.15833 82.45000];

% load static data
A = load('/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Processed/static_file_MERRA2.in');
lonlat = [A(:,3), A(:,2)];
lon = sort(unique(lonlat(:,1)));
lat = sort(unique(lonlat(:,2)));

% load fluxes
prec_names = dir('MERRA2*flx*.hdf');
rad_names = dir('MERRA2*rad*.hdf');
air_names = dir('MERRA2*slv*.hdf');

ndays = length(prec_names);
days_to_run = 1:ndays;
days_to_run = (ndays-10):ndays;

for k=[days_to_run] % do in parallel for individual MERRA-2 files?
    
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
    
    % Crop here, might cut down runtime significantly
    
    % interpolate to target resolution
    % this can use a lot of RAM depending on the target_resolution
    target_res = 1/16; % can change this if not enough RAM
    [prec_fine, target_lon, target_lat] = interpolate_merra2(prec, target_res, lon, lat); 
    temp_fine = interpolate_merra2(temp, target_res, lon, lat);
    ps_fine = interpolate_merra2(ps, target_res, lon, lat);
    vp_fine = interpolate_merra2(vp, target_res, lon, lat);
    wind_fine = interpolate_merra2(wind, target_res, lon, lat);
    swdn_fine = interpolate_merra2(swdn, target_res, lon, lat);
    lwdn_fine = interpolate_merra2(lwdn, target_res, lon, lat);
          
    %% crop to basin extent

    % create rectangle for cropping
    
    R1 = makerefmat(target_lon(1), target_lat(1), target_res, target_res);
    
    [ymin, xmin] = latlon2pix(R1, lat_range(1), lon_range(1));
    [ymax, xmax] = latlon2pix(R1, lat_range(2), lon_range(2));
    width = xmax - xmin;
    height = ymax - ymin;
    
    roundnumbers = 1;
    if roundnumbers
        rect = [floor(xmin), floor(ymin), ceil(width), ceil(height)];
    else
        rect = [xmin, ymin, width, height];
    end
    
    % crop each variable
    prec_crop = cell(24, 1);
    temp_crop = cell(24, 1);
    ps_crop = cell(24, 1);
    vp_crop = cell(24, 1);
    wind_crop = cell(24, 1);
    swdn_crop = cell(24, 1);
    lwdn_crop = cell(24, 1);
    for h=1:24
        prec_crop{h} = imcrop(prec_fine{h}, rect);
        temp_crop{h} = imcrop(temp_fine{h}, rect);
        ps_crop{h} = imcrop(ps_fine{h}, rect);
        vp_crop{h} = imcrop(vp_fine{h}, rect);
        wind_crop{h} = imcrop(wind_fine{h}, rect);
        swdn_crop{h} = imcrop(swdn_fine{h}, rect);
        lwdn_crop{h} = imcrop(lwdn_fine{h}, rect);
    end
%     clear prec_fine temp_fine ps_fine vp_fine wind_fine swdn_fine lwdn_fine % to conserve RAM
    
    out_lat = lat_range(1):target_res:lat_range(2);
    out_lon = lon_range(1):target_res:lon_range(2);
    nx = length(out_lon);
    ny = length(out_lat);
    R = makerefmat(out_lon(1), out_lat(1), target_res, target_res);

%     figure, imagesc(lon_range, lat_range, swdn_crop{12})
%     colorbar, xlabel('Lon'); ylabel('Lat');
%     set(gca, 'ydir', 'normal')

    %% Write VIC input files
    
    % go through the grid and print the 24 values for each grid cell to an
    % appropriately-named file
    for x=1:nx
        for y=1:ny          
            
            [row1, col1] = latlon2pix(R, out_lat(y), out_lon(x));

            prec_day = zeros(24, 1);
            temp_day = zeros(24, 1);
            ps_day = zeros(24, 1);
            vp_day = zeros(24, 1);
            wind_day = zeros(24, 1);
            swdn_day = zeros(24, 1);
            lwdn_day = zeros(24, 1);
            for h=1:24
                prec_day(h) = prec_crop{h}(row1, col1);
                temp_day(h) = temp_crop{h}(row1, col1);
                ps_day(h) = ps_crop{h}(row1, col1);
                vp_day(h) = vp_crop{h}(row1, col1);
                wind_day(h) = wind_crop{h}(row1, col1);
                swdn_day(h) = swdn_crop{h}(row1, col1);
                lwdn_day(h) = lwdn_crop{h}(row1, col1);
            end
            
            precision = '%3.5f';
            forcdir = '/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Forc'; % need to change this to a directory without spaces in the name...
            savename = fullfile(forcdir, ['Forcings_' num2str(out_lat(y), precision) '_' num2str(out_lon(x), precision) '.txt']);

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
            
            % write to file
            fID = fopen(savename, 'a');
            forcings_out = [temp_day, prec_day, ps_day, swdn_day, lwdn_day, vp_day, wind_day];
%             formatstring = '%14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n';
            formatstring = '%0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f\n';
            fprintf(fID, formatstring, forcings_out');
            fclose(fID);

        end
    end
    
    % Write dates to file
    % This is a check to make sure the VIC forcing files are in the
    % right order and aren't missing any data

    fID2 = fopen(fullfile(forcdir, 'forcing_dates.txt'), 'a');
    dates_out = [str2double(yr), str2double(mon), str2double(day)];
    formatstring = '%d %d %d\n';
    fprintf(fID2, formatstring, dates_out');
    fclose(fID2);
    
    disp(['Interpolated data for day ' num2str(k) ' of ' num2str(length(days_to_run))]) % progress tracker
    
end

%% Plot forcing time series for a particular grid cell

dat = load('/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/Shared_SWOTDA/MERRA2/Forc/Forcings_28.02500_66.40833.txt');

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










