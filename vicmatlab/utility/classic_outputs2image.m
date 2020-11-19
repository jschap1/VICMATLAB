% Classic (output) to Image
%
% Converts VIC classic mode outputs to image mode NetCDF format
% 8/18/2020 JRS
%
% INPUTS
% classic_outdir: location where classic mode VIC outputs are located
% image_outdir: save location for converted output file
% prefix: must include the underscore, if there is one. Example: "fluxes_"
% varnames: can specify variable names to convert, or leave blank to
% convert all variables
%
% Subsetting times has not been fully implemented yet. 

function ncoutname = classic_outputs2image(classic_outdir, prefix, image_outdir, varargin)

%% Initialization

numvarargs = length(varargin);
if numvarargs > 3
    error('The max number of optional arguments is 3')
end

% Option to convert all or just some of the outputs
optargs = {'all', NaN, NaN};
optargs(1:numvarargs) = varargin;
[varnames, start_date, end_date] = optargs{:};

if ~exist(image_outdir, 'dir')
    mkdir(image_outdir)
end

classic_names = dir(fullfile(classic_outdir, [prefix '*.txt']));
dat = dlmread(fullfile(classic_outdir, classic_names(1).name), '\t', 3, 0);
t = datetime(dat(:,1), dat(:,2), dat(:,3));

if isnan(start_date)    
    start_date = t(1);
    end_date = t(end);
else
    [~, i1] = ismember(t, start_date);
    [~, i2] = ismember(t, end_date);
    t = start_date:end_date;
end

if strcmp(varnames, 'all')
    varnames = get_vic_header(fullfile(classic_outdir, classic_names(1).name), 3);
end

ncoutname = fullfile(image_outdir, ...
    [prefix, num2str(year(start_date),'%02.f'), '-', num2str(month(start_date),'%02.f'), '-' num2str(day(start_date),'%02.f'), '.nc']);

%% Lon, lat, time, mask

[lon, lat] = get_coordinates_from_VIC_file(classic_outdir, prefix);
mask = xyz2grid(lon, lat, ones(length(lon), 1));
resolution = mode(abs(diff(lon)));

% if there is an issue with the grid, it's probably bc of this line. 
% Could change grid cell center vs. corner definition.
R = makerefmat(min(lon), min(lat), resolution, resolution); 

[lon1, lat1] = pixcenters(R, size(mask));

mask(isnan(mask)) = 0;
mask = single(mask);

nccreate(ncoutname,'lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon1)},...
          'Format','classic')      

nccreate(ncoutname,'lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat1)},...
          'Format','classic') 

nccreate(ncoutname,'time',...
    'Datatype','int32',...
    'Dimensions',{'time',length(t)},...
    'Format','classic')

nccreate(ncoutname,'mask',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon1),'lat',length(lat1)},...
          'Format','classic')
      
tt = (0:length(t)-1)'; % units are days since start_date
ncwrite(ncoutname,'time',tt);
ncwrite(ncoutname,'lon',lon1);
ncwrite(ncoutname,'lat',lat1);
ncwrite(ncoutname,'mask',mask');

ncwriteatt(ncoutname,...
    'time','calendar','proleptic_gregorian');
ncwriteatt(ncoutname,...
    'time','units',['days since ', num2str(year(start_date),'%02.f'), '-', num2str(month(start_date),'%02.f'), '-', num2str(day(start_date),'%02.f')]);

ncwriteatt(ncoutname,...
    'lon','standard_name','longitude');
ncwriteatt(ncoutname,...
    'lon','long_name','longitude of grid cell center');
ncwriteatt(ncoutname,...
    'lon','units','degrees_east');
ncwriteatt(ncoutname,...
    'lon','axis','X');

ncwriteatt(ncoutname,...
    'lat','standard_name','latitude');
ncwriteatt(ncoutname,...
    'lat','long_name','latitude of grid cell center');
ncwriteatt(ncoutname,...
    'lat','units','degrees_north');
ncwriteatt(ncoutname,...
    'lat','axis','Y');

ncwriteatt(ncoutname,...
    'mask','comment','NaN indicates cell is not active');
ncwriteatt(ncoutname,...
    'mask','long_name','fraction of grid cell that is activedomain mask');
ncwriteatt(ncoutname,...
    'mask','note','unitless');

%% Flux variables

nt = length(t);
nx = length(unique(lon));
ny = length(unique(lat));
for k=4:length(varnames)

    % Read in flux variable from classic mode outputs
    var = parload(classic_outdir, prefix, 0, k);
    
    % Make into a (masked) map to save
    output_map = NaN(ny, nx, nt); % this is the correct orientation
    for tt=1:nt
        output_map(:,:,tt) = xyz2grid(lon, lat, var(tt,:)');
    end
    
    % Removing the "OUT_" part of the name
    temp = strsplit(varnames{k}, '_');
    
    if length(temp) == 2
        varnames{k} = lower(temp{2});
    elseif length(temp) == 3
        varnames{k} = lower([temp{2} '_' temp{3}]);
    else
        varnames{k} = lower([temp{2} '_' temp{3}, '_', temp{4}]);
    end
    
    nccreate(ncoutname,varnames{k},...
        'Datatype','double',...
        'Dimensions',{'lon',length(lon1),'lat',length(lat1),'time',nt}, ...
        'Format','classic')

    ncwrite(ncoutname, varnames{k}, permute(output_map, [2,1,3]));

end

return



