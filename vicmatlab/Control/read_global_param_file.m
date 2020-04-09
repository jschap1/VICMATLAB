% Writes out the VIC global parameter file

function A = read_global_param_file(fname)    

fid = fopen(fname, 'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

return