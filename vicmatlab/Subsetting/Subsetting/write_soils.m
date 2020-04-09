% Write soils
% 
% Writes out the soil parameter file
% 
% INPUTS
% precision = number of decimal points for the coordinates/forcing file
% names
% soils = data for soil parameter file
% outname 
% outformat = livneh, 2l, or 3l

function write_soils(precision, soils, outname, outformat)

% Write out
fstring = ['%.' num2str(precision) 'f'];

if strcmp(outformat, 'livneh')
    % For the Livneh soil parameter file
    fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];
elseif strcmp(outformat, '2l')
    % For the 2-layer soil parameter file from HWSD
    fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %d %d %.3f %.3f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d\n'];
elseif strcmp(outformat, '3l')
    % This is used VIC-3L w no optional variables (54 columns in the soil parameter file).
    fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.4f %.4f %.2f %d %d %d %d %.2f\n'];
else
    error('Please specify a valid value of outformat')
end

fID = fopen(outname, 'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved to ' outname])


return