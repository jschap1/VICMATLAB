% Converts vegetation parameter files to geotiff files
%
% Updated for v1_3 of the VICGlobal parameters
%
% INPUTS
% Outputs from load_veg_parameters
% veglib = name of vegetation library file. Must be in numeric only format. Can leave blank.
%
% OUTPUTS
% geotiffs containing:
% -cover fraction (3)
% -rooting depths (3)
% -rooting fraction (3)
% -nveg (1)
% -landmask (1)

function [VEGPARAM] = convert_veg_parameters_v3(nvegtable, vegparamtable, latlontable, veglib, LC, outdir, plotflag, saveflag)

lat = latlontable(:,2);
lon = latlontable(:,3);

nlat = length(unique(lat));
nlon = length(unique(lon));

rasterSize = [nlat, nlon];
latlim = [min(lat), max(lat)];
lonlim = [min(lon), max(lon)];

R = georefcells(latlim,lonlim,rasterSize);

%% Create nveg map and landmask

ncells = length(lat);
mask_vect = nvegtable(:,2)>0;
nveg_map = xyz2grid(lon, lat, nvegtable(:,2));
vegmask = xyz2grid(lon, lat, mask_vect);

if plotflag
    figure, plotraster(lonlim, latlim, vegmask, 'Land mask from vegetation cover', 'Lon', 'Lat')
    figure, plotraster(lonlim, latlim, nveg_map, 'Number of vegetation classes', 'Lon', 'Lat')
end

if saveflag
    geotiffwrite(fullfile(outdir, 'nveg_map.tif'), nveg_map, R)
    geotiffwrite(fullfile(outdir, 'veg_mask.tif'), vegmask, R)
end

VEGPARAM.nveg = nveg_map;
VEGPARAM.landmask = vegmask;
VEGPARAM.R = R;
VEGPARAM.lat = lat;
VEGPARAM.lon = lon;
VEGPARAM.latlim = latlim;
VEGPARAM.lonlim = lonlim;

%% Create cover fraction maps for each vegetation type

classnames = fieldnames(vegparamtable);
nclasses = length(classnames);

for i=1:nclasses
    
     % start with i=16, water bodies
     % next, i=8, Evergreen_Needleleaf
     
     disp(['Making vegparam geotiffs for LC type: ', classnames{i}])
     
%     cellIDs = vegparamtable.(classnames{i})(:,1);
    vegfract = vegparamtable.(classnames{i})(:,2);
    rootdepth = vegparamtable.(classnames{i})(:,3:5);
    rootfract = vegparamtable.(classnames{i})(:,6:8);
    
    % get indices of cells where this vegetation type exists
    
    ind = find(vegfract > 0);
    vegfract_vect = zeros(ncells, 1);
    rootfract_arr = zeros(ncells, 3);
    rootdepth_arr = zeros(ncells, 3);
    vegfract_vect(ind) = vegfract(ind);
    rootfract_arr(ind, :) = rootfract(ind, :);
    rootdepth_arr(ind, :) = rootdepth(ind, :);

    vegfract_map = xyz2grid(lon, lat, vegfract_vect);
    
    rootfract1_map = xyz2grid(lon, lat, rootfract_arr(:,1));
    rootfract2_map = xyz2grid(lon, lat, rootfract_arr(:,2));
    rootfract3_map = xyz2grid(lon, lat, rootfract_arr(:,3));

    rootdepth1_map = xyz2grid(lon, lat, rootdepth_arr(:,1));
    rootdepth2_map = xyz2grid(lon, lat, rootdepth_arr(:,2));
    rootdepth3_map = xyz2grid(lon, lat, rootdepth_arr(:,3));
    
    if plotflag
        figure, plotraster(lonlim, latlim, vegfract_map, ['LC fraction ', classnames{i}], '', '')
        figure, plotraster(lonlim, latlim, rootdepth3_map, 'Root depth (root zone 3)', '', '')
        figure, plotraster(lonlim, latlim, rootfract3_map, 'Root fraction (root zone 3)', '', '')
    end    
    
    VEGPARAM.(classnames{i}).vegfract = vegfract_map;
    VEGPARAM.(classnames{i}).rootfract{1} = rootfract1_map;
    VEGPARAM.(classnames{i}).rootfract{2} = rootfract2_map;
    VEGPARAM.(classnames{i}).rootfract{3} = rootfract3_map;
    VEGPARAM.(classnames{i}).rootdepth{1} = rootdepth1_map;
    VEGPARAM.(classnames{i}).rootdepth{2} = rootdepth2_map;
    VEGPARAM.(classnames{i}).rootdepth{3} = rootdepth3_map;

if saveflag
    geotiffwrite(fullfile(outdir, [classnames{i}, '_vegfract.tif']), vegfract_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rf1.tif']), rootfract1_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rf2.tif']), rootfract2_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rf3.tif']), rootfract3_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rd1.tif']), rootdepth1_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rd2.tif']), rootdepth2_map, R)
    geotiffwrite(fullfile(outdir, [classnames{i}, '_rd3.tif']), rootdepth3_map, R)
end

 disp(['Finished processing rooting parameters for vegetation class ' num2str(i)])

end

%% % Plot each vegetation parameter from the vegetation library

% veglib = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_3/veglib_nh_nohead.txt';
% requires a special, modified version of the vegetation library without
% headers or comments -- makes it easier to load the vegetation library

monthnames = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

if ~isempty(veglib)
    
    vegetation_library = load(veglib);
    
     for i=1:nclasses
                
        disp(['Making veglib geotiffs for LC type: ', classnames{i}])
        
        % Initialize vegetation parameters for this land cover type
        overstory.(classnames{i}) = NaN(ncells, 1);
        Rarc.(classnames{i}) = NaN(ncells, 1);
        Rmin.(classnames{i}) = NaN(ncells, 1);
        for m=1:12
            LAI.(monthnames{m}).(classnames{i}) = NaN(ncells, 1);
            Fcan.(monthnames{m}).(classnames{i}) = NaN(ncells, 1);
            albedo.(monthnames{m}).(classnames{i}) = NaN(ncells, 1);
        end
        roughness.(classnames{i}) = NaN(ncells, 1);
        displacement.(classnames{i}) = NaN(ncells, 1);
        wind_h.(classnames{i}) = NaN(ncells, 1);
        RGL.(classnames{i}) = NaN(ncells, 1);
        
        % get index for cells of this land cover type with cover fraction > 0
        ind1 = vegparamtable.(classnames{i})(:,2) > 0;
        
        % get the values from the vegetation library
        % LC.classnumber(i) converts from alphabetical order (i=1:17) to
        % the order of the IGBP indexes in the vegetation library file
        veglibvals.overstory = vegetation_library(i, 2);
        veglibvals.Rarc = vegetation_library(i, 3);
        veglibvals.Rmin = vegetation_library(i, 4);
        for m=1:12
            veglibvals.LAI.(monthnames{m}) = vegetation_library(i, 4+m);
            veglibvals.Fcan.(monthnames{m}) = vegetation_library(i, 16+m);
            veglibvals.albedo.(monthnames{m}) = vegetation_library(i, 28+m);
        end
        % any month will do for roughness and displacement height, since
        % they are assumed temporally invariant
        veglibvals.roughness = vegetation_library(i, 41);
        veglibvals.displacement = vegetation_library(i, 53);
        veglibvals.wind_h = vegetation_library(i, 65);
        veglibvals.RGL = vegetation_library(i, 66);
        % later, I want to use IGBP numbers for the vegetation classes,
        % instead of alphabetical order, but this will do for now
        
        % set values of the overstory variable for this vegetation type
        overstory.(classnames{i})(ind1) = veglibvals.overstory;
        Rarc.(classnames{i})(ind1) = veglibvals.Rarc;
        Rmin.(classnames{i})(ind1) = veglibvals.Rmin;
        for m=1:12
            LAI.(monthnames{m}).(classnames{i})(ind1) = veglibvals.LAI.(monthnames{m});
            Fcan.(monthnames{m}).(classnames{i})(ind1) = veglibvals.Fcan.(monthnames{m});
            albedo.(monthnames{m}).(classnames{i})(ind1) = veglibvals.albedo.(monthnames{m});
        end
        roughness.(classnames{i})(ind1) = veglibvals.roughness;
        displacement.(classnames{i})(ind1) = veglibvals.displacement;
        wind_h.(classnames{i})(ind1) = veglibvals.wind_h;
        RGL.(classnames{i})(ind1) = veglibvals.RGL;
        
        % reshape the vegetation parameter vectors into gridded maps
        overstory.map.(classnames{i}) = xyz2grid(lon, lat, overstory.(classnames{i}));
        Rarc.map.(classnames{i}) = xyz2grid(lon, lat, Rarc.(classnames{i}));
        Rmin.map.(classnames{i}) = xyz2grid(lon, lat, Rmin.(classnames{i}));
        for m=1:12
            LAI.(monthnames{m}).map.(classnames{i}) = xyz2grid(lon, lat, LAI.(monthnames{m}).(classnames{i}));
            Fcan.(monthnames{m}).map.(classnames{i}) = xyz2grid(lon, lat, Fcan.(monthnames{m}).(classnames{i}));
            albedo.(monthnames{m}).map.(classnames{i}) = xyz2grid(lon, lat, albedo.(monthnames{m}).(classnames{i}));
        end        
        roughness.map.(classnames{i}) = xyz2grid(lon, lat, roughness.(classnames{i}));
        displacement.map.(classnames{i}) = xyz2grid(lon, lat, displacement.(classnames{i}));
        wind_h.map.(classnames{i}) = xyz2grid(lon, lat, wind_h.(classnames{i}));
        RGL.map.(classnames{i}) = xyz2grid(lon, lat, RGL.(classnames{i}));
        
        if plotflag
            figure, plotraster(lonlim, latlim, overstory.map.(classnames{i}), ['Overstory ' classnames{i}], '', '')
            figure, plotraster(lonlim, latlim, LAI.Jan.map.(classnames{i}), ['January LAI ' classnames{i}], '', '')
            figure, plotraster(lonlim, latlim, LAI.Jul.map.(classnames{i}), ['July LAI ' classnames{i}], '', '')
            figure, plotraster(lonlim, latlim, RGL.map.(classnames{i}), ['RGL ' classnames{i}], '', '')
        end
       
        if saveflag
            
            % overstory
            outname = fullfile(outdir, ['overstory_', classnames{i}, '.tif']);
            geotiffwrite(outname, overstory.map.(classnames{i}), R)
            
            % architectural resistance
            outname = fullfile(outdir, ['Rarc_', classnames{i}, '.tif']);
            geotiffwrite(outname, Rarc.map.(classnames{i}), R)
            
            % minimum radiation for transpiration
            outname = fullfile(outdir, ['Rmin_', classnames{i}, '.tif']);
            geotiffwrite(outname, Rmin.map.(classnames{i}), R)
                
            for m=1:12
                
                % LAI
                outname = fullfile(outdir, ['LAI_', classnames{i} '_' monthnames{m}, '.tif']);
                geotiffwrite(outname, LAI.(monthnames{m}).map.(classnames{i}), R)
                
                % canopy fraction
                outname = fullfile(outdir, ['Fcan_', classnames{i} '_' monthnames{m}, '.tif']);
                geotiffwrite(outname, Fcan.(monthnames{m}).map.(classnames{i}), R)                
                
                % albedo
                outname = fullfile(outdir, ['albedo_', classnames{i} '_' monthnames{m}, '.tif']);
                geotiffwrite(outname, albedo.(monthnames{m}).map.(classnames{i}), R)
                
            end
            
            % roughness length
            outname = fullfile(outdir, ['roughness_', classnames{i}, '.tif']);
            geotiffwrite(outname, roughness.map.(classnames{i}), R)
            
            % displacement height
            outname = fullfile(outdir, ['displacement_', classnames{i}, '.tif']);
            geotiffwrite(outname, displacement.map.(classnames{i}), R)
            
            % wind measurement height
            outname = fullfile(outdir, ['wind_h', classnames{i}, '.tif']);
            geotiffwrite(outname, wind_h.map.(classnames{i}), R)
            
            % RGL
            outname = fullfile(outdir, ['RGL_', classnames{i}, '.tif']);
            geotiffwrite(outname, RGL.map.(classnames{i}), R)            

        end
        
        disp(['Finished processing vegetation library parameters for vegetation class ' num2str(i)])
        
    end

    % Save output
    disp('Saving processed vegetation parameters');
    
    for i=1:nclasses
    
        VEGPARAM.(classnames{i}).overstory = overstory.(classnames{i});
        VEGPARAM.(classnames{i}).Rarc = Rarc.(classnames{i});
        VEGPARAM.(classnames{i}).Rmin = Rmin.(classnames{i});
        
        for m=1:12
            VEGPARAM.(classnames{i}).LAI.(monthnames{m}) = LAI.(monthnames{m}).(classnames{i});
            VEGPARAM.(classnames{i}).Fcan.(monthnames{m}) = Fcan.(monthnames{m}).(classnames{i});
            VEGPARAM.(classnames{i}).albedo.(monthnames{m}) = albedo.(monthnames{m}).(classnames{i});
        end       
        
        VEGPARAM.(classnames{i}).roughness = roughness.(classnames{i});  
        VEGPARAM.(classnames{i}).displacement = displacement.(classnames{i});
        VEGPARAM.(classnames{i}).wind_h = wind_h.(classnames{i});
        VEGPARAM.(classnames{i}).RGL = RGL.(classnames{i});
        
    end
    
    save(fullfile(outdir, 'VEGPARAM.mat'), 'VEGPARAM', '-v7.3')
    disp(['Saved processed vegetation parameters to ' fullfile(outdir, 'VEGPARAM.mat')])
    
end

return

