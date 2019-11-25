% Plot one grid cell
%
% Plots the time series of a VIC output at a user-specified location
%
% Inputs
% lat, lon = coordinates where you would like to plot the output variable
% info = structure containing metadata from the VIC run. See get_vic_run_metadata()
% ebwbflag = 'eb' or 'wb'
% variable_name = name of variable you'd like to plot. Must be an exact match.
%
% Outputs
% f = figure handle

function f = plot_one_grid_cell(lat, lon, info, ebwbflag, variable_name)

T = delaunayn([info.lat, info.lon]);
k = dsearchn([info.lat, info.lon],T,[lat lon]);

if strcmp(ebwbflag, 'wb')
    
    wbname = fullfile(info.wb_out_dir, ['wb_' num2str(info.lat(k),8) '_' num2str(info.lon(k),8) '.txt']);
    wb_out = dlmread(wbname, '\t', info.headerlines, 0);
    wb_out = array2table(wb_out);
    wb_out.Properties.VariableNames = info.wbvars; 
    f = jsplot(info.time, wb_out.(variable_name), variable_name, 'Time', variable_name, 14);
    title(variable_name, 'Interpreter', 'none');
    
elseif strcmp(ebwbflag, 'eb')
    
    ebname = fullfile(info.eb_out_dir, ['eb_' num2str(info.lat(k),8) '_' num2str(info.lon(k),8) '.txt']);
    eb_out = dlmread(ebname, '\t', info.headerlines, 0);
    eb_out = array2table(eb_out);
    eb_out.Properties.VariableNames = info.ebvars;
    f = jsplot(info.time, eb_out.(variable_name), variable_name, 'Time', variable_name, 14);
    title(variable_name, 'Interpreter', 'none');
    
else
    disp('Enter a valid value for ebwbflag')
end

return