% Creates input files for the Lohmann routing model
%
% Written 10/23/2019 JRS
% Updated (made function) 8/12/2020 JRS
% Could be improved re: reading in VIC output files; there is some relevant
% code elsewhere in VICMATLAB
%
% Also useful:
% Remove extension from VIC output files
% for file in fluxes*.txt; do mv "$file" "${file%.txt}"; done

function convert_routing_inputs(vic_out_dir, rout_in_dir)

fnames = dir([vic_out_dir '/wb_*']);

ncells = length(fnames);

if ~exist(rout_in_dir, 'dir')
    mkdir(rout_in_dir)
    disp(['Created directory: ' rout_in_dir])
end

for k=1:ncells
    
    dat = dlmread(fullfile(vic_out_dir, fnames(k).name), '\t', 3, 0);
    nt = size(dat,1);
%     A = [dat(:,1:3), zeros(nt,1), zeros(nt,1), dat(:,6), dat(:,7)];
    
    A = [dat(:,1:3), zeros(nt,1), zeros(nt,1), dat(:,5), dat(:,6)];
    outname = ['fluxes_' fnames(k).name(4:end)]; % if name is wb_
%     outname = ['fluxes_' fnames(k).name(10:end)]; % if name is fluxes_ or wb_daily
    dlmwrite(fullfile(rout_in_dir, outname), A, '\t')
    
end

return