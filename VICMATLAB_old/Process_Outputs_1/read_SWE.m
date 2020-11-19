function [swe, t, ds] = read_SWE(indir, ncells, nt)

ds = tabularTextDatastore(indir);
t = tall(ds);

ds.ReadSize = 'file';

reset(ds)
swe = zeros(ncells, nt);

for k = 1:ncells
    T = read(ds);
    swe(k, :) = table2array(T(:,27));
end

return