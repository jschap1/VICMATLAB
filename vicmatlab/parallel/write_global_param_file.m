% Writes out the VIC global parameter file

function write_global_param_file(A, outname)    

fid = fopen(outname, 'w');
for i=1:numel(A)
    if A{i+1} == -1 % MATLAB may read the final line as '-1'
        fprintf(fid,'%s\n', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);

return