% Save Forcings
% 
% 4/2/2020 JRS
% A function to save the outputs of subset_forcings

function save_forcings(metlat, lat_ind, metlon, lon_ind, force_out, data_cum, grid_decimal)

ncells = length(lat_ind);
fstring = ['%.' num2str(grid_decimal) 'f'];
for k=1:ncells     
    filename = ['data_' num2str(metlat(lat_ind(k)),fstring) '_' num2str(metlon(lon_ind(k)),fstring)];
    dlmwrite(fullfile(force_out, filename), data_cum(:,:,k), ' ')
end
disp(['Forcing data saved to ' force_out])

end