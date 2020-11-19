% Reassign soil parameters
% 
% 8/16/2020 JRS
% Reassigns soil parameters from one soil parameter file (e.g. L2013) to
% another (e.g. VICGlobal).
%
% INPUTS
% spf1 = name of first soil parameter file
% spf2 = name of second soil parameter file
% cols = columns of first soil parameter file to transfer to second
%
% Note: code requires columns of first and second spfs to match

function soils_VG = reassign_soil_parameters(spf1, spf2, cols)

soils_L13 = load(spf1);
soils_VG = load(spf2);

lat_vg = soils_VG(:,3);
lon_vg = soils_VG(:,4);

lat_L13 = soils_L13(:,3);
lon_L13 = soils_L13(:,4);

T = delaunayn([lat_vg, lon_vg]);
ncells = length(lon_vg);
disp(['There are ' num2str(ncells) ' grid cells'])
nvars = length(cols);
for j=1:ncells
    k = dsearchn([lat_L13, lon_L13],T,[lat_vg(j) lon_vg(j)]);
%     soils_VG(j,3:4) == soils_L13(k,3:4); % lon/lat
    for i = 1:nvars
        soils_VG(j,cols(i)) = soils_L13(k,cols(i));
    end
end

outformat = '3l';
precision = 5;
temp = strsplit(spf2, '.txt');
outname = [temp{1} '_reassigned.txt'];
write_soils(precision, soils_VG, outname, outformat)

return