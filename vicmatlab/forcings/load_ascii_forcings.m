% Load ASCII Forcing Data
%
% Loads ASCII forcing data (e.g. for plotting)
%
% Dependencies: GetCoords.m
% Updated 4/2/2020 JRS
%
% INPUTS
% forcingpath = location of ASCII forcing files to plot
% prefix
% precision
% varnames
%
% OUTPUTS
% forc = structure array with information about the forcing data
%
% EXAMPLES
% forcingpath = '/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Forcings/Disagg_Forc/'; 
% forcingpath = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain2/aligned_forcings';
% forcingpath = '/Volumes/HD4/SWOTDA/Data/Tuolumne/forc_ascii';
% forcingpath = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain2/2018-2018_forc_lasttest';
%
% Sample arguments
% forcingpath = './data/forc_2009-2011';
% precision = 5; 
% varnames = {'PRECIP','TMIN','TMAX','WIND'};
% prefix = 'data_';

function forc = load_ascii_forcings(forcingpath, prefix, precision, varnames)

%% Load forcing data

forcenames = dir(fullfile(forcingpath, [prefix '*']));
ncells = length(forcenames);
tmpforc = dlmread(fullfile(forcingpath, forcenames(1).name)); 

nsteps = size(tmpforc,1);
nvars = size(tmpforc,2);

% Adapted from LoadVICResults
gridcells = cell(ncells, 1);
for k=1:ncells
    tmpstring = forcenames(k).name;
%     tmpstring = strrep(tmpstring,'-',''); % remove some characters bc Matlab cannot handle them
    tmpstring = strrep(tmpstring,'.','_');
    gridcells{k} = tmpstring;
end

[lat, lon] = GetCoords(gridcells, precision, prefix);

%% Read in data for all grid cells and times (OK for a relatively small domain)

% This is definitely not going to work for a large domain
% Going to have to break it up
% Check if the domain is "large"

if ncells > 500
    disp('Relatively large study area')
    disp('Use make_netcdf_forcing()')
    error('Study area is too big. It will fill up RAM')
%     disp('Loading daily average forcing data')
%     
%     % Using a matfile to deal with the memory issue
%     matOjb = matfile('./precipitation.mat');
%     
%     ncells = 10;
% %     precip = NaN(nsteps, ncells);
% 
%     for k=1:ncells
%         forc_current = dlmread(fullfile(forcingpath, forcenames(k).name));
%         
%         precip = forc_current(:,1);
%                 
%         if k==1
%             save('./precip.m', 'precip', '-v7.3')
%             m = matfile('./precip.m', 'Writable', true);
%         else
%             m.precip(:,k) = forc_current;
%         end
%         
%     end
    
else
    disp('Relatively small study area')
    disp('Loading all forcing data')
    FORC = NaN(nsteps, nvars, ncells);
    disp(['There are ' num2str(nvars) ' variables in the forcing files']);
    for k=1:ncells
        FORC(:,:,k) = dlmread(fullfile(forcingpath, forcenames(k).name)); 
    end    
end


%% Put forcing data into a nice structure

forc = struct();
forc.names = varnames;
for i=1:length(varnames)
    forc.(varnames{i}) = squeeze(FORC(:,i,:));
end
forc.dimensions = {'timestep','gridcell'};
forc.lon = lon;
forc.lat = lat;

end
