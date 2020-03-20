% Removes .txt extension from VIC5 flux output files
% This is necessary to get the Lohmann routing model to work - the legacy
% model written in Fortran and included with VIC4.2
%
% This takes a long time when there are a lot of VIC output files.
% Like, weeks.
%
% But it can be done in like a second using shell scripting:
% for file in fluxes*.txt; do mv "$file" "${file%.txt}"; done

orignames = dir('fluxes*');
ncells = length(orignames);

for k=1:ncells
    ext_ind = strfind(orignames(k).name,'.txt');
    newname = orignames(k).name(1:ext_ind-1);
    movefile(orignames(k).name,newname);
    delete(orignames(k).name);
    k
end

orignames = dir('Rout*');
ncells = length(orignames);

for k=1:ncells
    ext_ind = strfind(orignames(k).name,'R_Out');
    newname = orignames(k).name(13:end);
    movefile(orignames(k).name,newname);
    delete(orignames(k).name);
    k
end

%% Remove headers

% VIC 5 Classic output files come with three lines of header rows. In order
% to run the VIC routing model with VIC 5 outputs, these header rows need
% to be removed.

in_dir = './Outputs/VIC_IRB/WB_corrected_veg_Raw';
out_dir = './Outputs/VIC_IRB/WB_corrected_veg/v4';

orignames = dir(fullfile(in_dir, 'fluxes*'));
ncells = length(orignames);
header_rows = 3;
for k=1:ncells
    data = dlmread(fullfile(in_dir, orignames(k).name), '\t', header_rows, 0);
    dlmwrite(fullfile(out_dir, orignames(k).name), data, ' ');
    if mod(k, 1000)==0
        disp(k)
    end
end



