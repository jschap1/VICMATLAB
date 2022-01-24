% Upscale vegetation parameter file
% 
% Reads in the vegetation parameter file
% Averages parameter values, aggregating them to a higher (or lower) resolution
%
% Soils must be at same resolution as orig_vegpar

function upscale_vegpars(vegparfile, soils_old, soils_new, newres, oldres)

lon = soils_old(:,4);
lat = soils_old(:,3);
elevmap = xyz2grid(lon, lat, soils_old(:,22));
figure
plotraster(lon, lat, elevmap,'Elevation (m)')

% Read in vegetation parameter file
ncells = size(soils_old,1);
savename = './Data/CONUS_veg.mat';

if exist(savename, 'file') == 0
    [nvegtable, vegparamtable, latlontable, LC] = read_vegparam_w_LAI(vegparfile, soils_old, ncells, savename);
else
    load(savename);
end

% f = newres/oldres; % upscaling (or downscaling) factor
[n, nvars] = size(soils_old); % n = original number of grid cells

% [m1, n1] = size(elevmap); % original dimensions

% new dimensions
% m2 = m1/f;
% n2 = floor(n1/f);

% Get (any) field to upscale to initialize the method
LC = fieldnames(vegparamtable);
invar = vegparamtable.Closed_Shrublands(:,2); % cover fraction
invarmap = xyz2grid(latlontable(:,2), latlontable(:,3), invar)';

figure
plotraster(latlontable(:,2), latlontable(:,3), invarmap, 'Cover fraction')

[outvar, newlons, newlats] = upscale_raster(invarmap, lon, lat, newres, oldres, 'nearest');
% outvar = interp2(oldlats, oldlons, invarmap', newlats, newlons, 'nearest')';

figure
plotraster(newlons, newlats, outvar, 'Cover fraction')

%% Cover fraction

n = size(soils_new,1); % number of grid cells in output
cellID = soils_new(:,2);

% Cover fraction
cover_fraction = struct();
for k=1:length(LC) % first variable is cellID, so start at k=2
    
    invar = vegparamtable.(LC{k})(:,2); % cover fraction
    invarmap = xyz2grid(latlontable(:,2), latlontable(:,3), invar)';
    
    cover_fraction.(LC{k}) = upscale_raster(invarmap, lon, lat, newres, oldres, 'nearest');
%     cover_fraction.(LC{k}) = interp2(oldlats, oldlons, invarmap', newlats, newlons, 'nearest')';
    
end

% Prepare for masking (could be done more efficiently using soils input)
runcell = xyz2grid(lon, lat, soils_old(:,1));
runcell2 = upscale_raster(runcell, lon, lat, newres, oldres, 'linear');
% runcell2 = runcell2';
% runcell2(isnan(runcell2))=0;
runcell3 = runcell2(:)==0;

gridmask = xyz2grid(soils_new(:,4), soils_new(:,3), soils_new(:,1));
figure
plotraster(newlons, newlats, gridmask, 'mask');

% Calculate nveg
nveg = zeros(size(outvar));
for k=1:length(LC)
    tmp = cover_fraction.(LC{k});
    tmp(isnan(tmp)) = 0;
    nveg = nveg + logical(tmp);
end
nveg = nveg(:);
nveg(runcell3) = [];

% Calculate vegtypes
vgn = zeros(n,length(LC));
for k=1:length(LC)
    tmp = cover_fraction.(LC{k});
%     tmp(isnan(tmp)) = 0;
    tmp = tmp(:);
    tmp(runcell3) = [];
    vgn(:,k) = tmp > 0;
end
vgn2 = cell(n,1);
for kk=1:n
    vgn2{kk} = find(vgn(kk,:))';
end

%% Root fraction and depth

% header_names = {'cover_fraction','rootfract1','rootfract2','rootfract3', 'rootdepth1','rootdepth2', 'rootdepth3', ...
%     'LAI1','LAI2','LAI3','LAI4','LAI5','LAI6','LAI7','LAI8','LAI9','LAI10','LAI11','LAI12'};

rd1_str = struct();
rd2_str = struct();
rd3_str = struct();
rf1_str = struct();
rf2_str = struct();
rf3_str = struct();

for k=1:length(LC) % first variable is cellID, so start at k=2
    
    rd1 = vegparamtable.(LC{k})(:,3);
    rf1 = vegparamtable.(LC{k})(:,4);
    rd2 = vegparamtable.(LC{k})(:,5);
    rf2 = vegparamtable.(LC{k})(:,6);
    rd3 = vegparamtable.(LC{k})(:,7);
    rf3 = vegparamtable.(LC{k})(:,8);
    
    rd1_map = xyz2grid(latlontable(:,2), latlontable(:,3), rd1)';
    rd2_map = xyz2grid(latlontable(:,2), latlontable(:,3), rd2)';
    rd3_map = xyz2grid(latlontable(:,2), latlontable(:,3), rd3)';
    rf1_map = xyz2grid(latlontable(:,2), latlontable(:,3), rf1)';
    rf2_map = xyz2grid(latlontable(:,2), latlontable(:,3), rf2)';
    rf3_map = xyz2grid(latlontable(:,2), latlontable(:,3), rf3)';
    
    rf1_str.(LC{k}) = upscale_raster(rf1_map, lon, lat, newres, oldres, 'nearest');
    rf2_str.(LC{k}) = upscale_raster(rf2_map, lon, lat, newres, oldres, 'nearest');
    rf3_str.(LC{k}) = upscale_raster(rf3_map, lon, lat, newres, oldres, 'nearest');
    
    rd1_str.(LC{k}) = upscale_raster(rd1_map, lon, lat, newres, oldres, 'nearest');
    rd2_str.(LC{k}) = upscale_raster(rd2_map, lon, lat, newres, oldres, 'nearest');
    rd3_str.(LC{k}) = upscale_raster(rd3_map, lon, lat, newres, oldres, 'nearest');
    
%     rf1_str.(LC{k}) = interp2(oldlats, oldlons, rf1_map', newlats, newlons, 'nearest')';
%     rf2_str.(LC{k}) = interp2(oldlats, oldlons, rf2_map', newlats, newlons, 'nearest')';
%     rf3_str.(LC{k}) = interp2(oldlats, oldlons, rf3_map', newlats, newlons, 'nearest')';
    
%     rd1_str.(LC{k}) = interp2(oldlats, oldlons, rd1_map', newlats, newlons, 'nearest')';
%     rd2_str.(LC{k}) = interp2(oldlats, oldlons, rd2_map', newlats, newlons, 'nearest')';
%     rd3_str.(LC{k}) = interp2(oldlats, oldlons, rd3_map', newlats, newlons, 'nearest')';
    
    % Check
%     rf = rf1_str.Evergreen_Needleleaf + rf2_str.Evergreen_Needleleaf + rf3_str.Evergreen_Needleleaf;
%     rf_orig = rf1_map + rf2_map + rf3_map;
%     figure,plotraster(lon, lat, rf_orig, '')
%     figure,plotraster(lon, lat, rf, '')
    
end

root_depth = cell(length(LC), 1);
root_fraction = cell(length(LC), 1);
for k=1:length(LC)
    
    tmp1 = rf1_str.(LC{k})(:);
    tmp1(runcell3==1) = [];
    tmp2 = rf2_str.(LC{k})(:);
    tmp2(runcell3==1) = [];
    tmp3 = rf3_str.(LC{k})(:);
    tmp3(runcell3==1) = [];
    root_fraction{k} = [tmp1, tmp2, tmp3];
    
    tmp1 = rd1_str.(LC{k})(:);
    tmp1(runcell3==1) = [];
    tmp2 = rd2_str.(LC{k})(:);
    tmp2(runcell3==1) = [];
    tmp3 = rd3_str.(LC{k})(:);
    tmp3(runcell3==1) = [];
    root_depth{k} = [tmp1, tmp2, tmp3];
    
end
% goal is to populate the output vegpar file

%% LAI

LAI_str = struct();

for k=1:length(LC) % first variable is cellID, so start at k=2
    
    LAI1 = vegparamtable.(LC{k})(:,9);
    LAI2 = vegparamtable.(LC{k})(:,10);
    LAI3 = vegparamtable.(LC{k})(:,11);
    LAI4 = vegparamtable.(LC{k})(:,12);
    LAI5 = vegparamtable.(LC{k})(:,13);
    LAI6 = vegparamtable.(LC{k})(:,14);
    LAI7 = vegparamtable.(LC{k})(:,15);
    LAI8 = vegparamtable.(LC{k})(:,16);
    LAI9 = vegparamtable.(LC{k})(:,17);
    LAI10 = vegparamtable.(LC{k})(:,18);
    LAI11 = vegparamtable.(LC{k})(:,19);
    LAI12 = vegparamtable.(LC{k})(:,20);
    
    LAI1_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI1)';
    LAI2_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI2)';
    LAI3_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI3)';
    LAI4_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI4)';
    LAI5_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI5)';
    LAI6_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI6)';
    LAI7_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI7)';
    LAI8_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI8)';
    LAI9_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI9)';
    LAI10_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI10)';
    LAI11_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI11)';
    LAI12_map = xyz2grid(latlontable(:,2), latlontable(:,3), LAI12)';

    for mm=1:12
%         LAI_str.(LC{k}).(['LAI' num2str(mm)]) = ...
%             interp2(oldlats, oldlons, eval(['transpose(LAI' num2str(mm) '_map)']), newlats, newlons, 'nearest')';
%         
            LAI_str.(LC{k}).(['LAI' num2str(mm)]) = ...
                upscale_raster(eval(['LAI' num2str(mm) '_map']), lon, lat, newres, oldres, 'nearest');
    
    end
    disp(k)
    
end

LAI = cell(length(LC), 1);
for k=1:length(LC)
    
    LAIvals = zeros(n,mm);
    for mm=1:12
        tmp = LAI_str.(LC{k}).(['LAI' num2str(mm)])(:);
        tmp(runcell3==1) = [];
        LAIvals(:,mm) = tmp;
    end

    LAI{k} = LAIvals;
    
end

%% Write out the upscaled vegetation parameter file

outname = './Data/CONUS_vegpars_1_4.txt';
fID = fopen(outname, 'w');

for k=1:n

    fmt = '%d %d\n';
    fprintf(fID, fmt, [cellID(k), nveg(k)]);
  
    % For each vegetation class, write parameters, then (optionally) monthly LAI

    for ii=1:nveg(k)
        current_vgn = vgn2{k}(ii);
        
        % Clunky to do this in the loop, but it works
        tmp = cover_fraction.(LC{current_vgn})(:);
        tmp(runcell3==1) = [];
        current_cf = tmp(k);
        
        current_rd = root_depth{current_vgn}(k,:);
        current_rf = root_fraction{current_vgn}(k,:);
        current_LAI = LAI{current_vgn}(k,:);
        
        fmt = '%d %4.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n';
        
%         vegpars = [current_vgn current_cf, current_rd, current_rf];
        vegpars = [current_vgn current_cf, current_rd(1), current_rf(1), ...
            current_rd(2), current_rf(2), current_rd(3), current_rf(3)];
        
        fprintf(fID, fmt, vegpars);
                
        fmt = '%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n';
       
        current_LAI(current_LAI<=0.01) = 0.01; % no zero LAI allowed
        fprintf(fID, fmt, current_LAI);
        
    end
               
end

fclose(fID);


return

