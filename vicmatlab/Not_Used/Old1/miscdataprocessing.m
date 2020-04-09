%%

% One variable, 16 grid cells, one time

p = 1;
t = 1; 

prec = NaN(ncells,1);
evap = NaN(ncells,1);
for k=1:ncells
    
    p=1;
    prec(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
    
    p=2;
    evap(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
    
end

FLUXES.maps.lat = lat;
FLUXES.maps.lon = lon;

(varnames{1})


%%

% Twenty variables, 16 grid cells, one time
% This is really all that is necessary.

for p = 1:length(varnames)
    
    FLUXES.maps.(varnames{p}) = NaN(ncells,1);
    
    for k=1:ncells

        FLUXES.maps.(varnames{p})(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
        % Currently only takes the first soil moisture layer. Modify to
        % split up each soil moisture layer into its own variable.
        
    end
    
end


%%

% Twenty variables, 16 grid cells, ntimesteps

for t = 1:ntimesteps
    
    for p = 1:length(varnames)

        FLUXES.maps.(varnames{p}) = NaN(ncells,1);

        for k=1:ncells

            FLUXES.maps.(varnames{p})(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
            % Currently only takes the first soil moisture layer. Modify to
            % split up each soil moisture layer into its own variable.

        end

    end

end



%%





% Dynamic fieldnames:

handles.button1 = 10;
handles.button2 = 20;
handles.button3 = 30;

for i=1:n
    control_name = ['button' num2str(i)];
    disp(handles.(control_name)); % this is a good way to use dynamic field names
end
%%
varnames = FLUXES.ts.gridcell_fluxes_48_1875_120_6875.Properties.VariableNames;

all_cells = cell(20,1);
for k = 1:ncells
    for p = 1:length(varnames)
        %FLUXES.maps(1).(varnames{p}) = FLUXES.ts.(cellnames{k}).(varnames{p})(1,:);
        
        
        all_cells{k} = FLUXES.ts.(cellnames{k}).(varnames{p})(1,:);
        
        
        %FLUXES.maps(1).(varnames{p}) = 
        
        
    end
end

%%
%cmd = ['FLUXES.maps(t).' varname{p} ' = FLUXES.ts.' cellnames{k} '.' varnames{p} '(t);'];

cellnames = fieldnames(FLUXES.ts);

flux_all_cells = NaN(ncells,1);
for k=1:ncells
    cmd = ['flux_all_cells(k) = FLUXES.ts.' cellnames{k} '.prec(1)'];
    eval(cmd);
end
FLUXES.maps(1).prec = flux_all_cells;

%FLUXES.ts.gridcell_fluxes_48_1875_120_6875.prec(1);

%%

lat
lon
ncells
cellnames
ntimesteps = 18993;
varnames

fluxes.ts = NaN(ncells, ntimesteps);

['FLUXES.ts.' cellnames{k} '.prec;']

for t=1:ntimesteps
tic
for t=1:ntimesteps
    for k=1:ncells
        cmd = ['FLUXES.maps(t) = FLUXES.ts.' cellnames{k} '.' varnames{p} '(1);']
        eval(cmd);
    end
end

for p=1:length(varnames);
    for k=1:ncells
        cmd = ['fluxes.ts(k, p) = FLUXES.ts.' cellnames{k} '.' varnames{p} '(1);'];
        eval(cmd);
    end
end


toc

% The result is a table with lat, lon, and corresponding flux values. 
% This can be plotted as a map relatively easily, and it can also be
% converted to a raster and georeferenced.
xyz = [lat, lon, fluxes.ts];
FLUXES.maps.prec = table(xyz);