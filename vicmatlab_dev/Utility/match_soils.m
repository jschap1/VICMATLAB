soils_VG = load('/Volumes/HD4/SWOTDA/Data/UMRB/soils_umrb_vg.txt');
soils_L15 = load('/Volumes/HD4/SWOTDA/Data/UMRB/soils_umrb_livneh.txt');

masklat = soils_L15(:,3);
masklon = soils_L15(:,4);

soils_VG(:,1) = 0;
T = delaunayn([soils_VG(:,3), soils_VG(:,4)]);
for i = 1:size(masklat,1)
    k = dsearchn([soils_VG(:,3), soils_VG(:,4)],T,[masklat(i) masklon(i)]);
    soils_VG(k,1) = 1;
    if mod(i,1e3) == 0
        disp(i)
    end
end

sum(soils_VG(:,1))

write_soils(5, soils_VG, '/Volumes/HD4/SWOTDA/Data/UMRB/soils_umrb_vg_2.txt', '3l')