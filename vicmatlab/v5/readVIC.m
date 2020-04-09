% Read VIC output files

function var_out_all = readVIC(info, fnames, varnames)

formatspec = [repmat('%f ', 1, length(varnames)-1) '%f\n'];
% wb_out_all = NaN(info.nt, length(info.wb), ncells); % time, variables, cells
% sizeA = [info.nt, length(info.wb)]; % size of the array being read in by fscanf

var_out_all = NaN(info.nt, length(varnames), info.ncells);

% tic
% This could be adapted to read in a particular variable
% Best practice is probably to save the outputs to disk as it goes
% so it doesn't need to use RAM
% Could use the Matlab big data functions for this
for k=1:info.ncells

    fID = fopen(fullfile(info.wb_out_dir, fnames(k).name));
    % skip the headerlines
    for kk=1:info.headerlines
        fgetl(fID);
    end    
    var_out_all(:,:,k) = cell2mat(textscan(fID, formatspec));
    fclose(fID);
    % textscan method takes 16.2 seconds for 1000 files
    
    % dlmread slow (takes 27.92 seconds for 1000 files)
    
%     wb_out_all(:,:,k) = dlmread(fullfile(vic_out_dir, wbnames(k).name), '\t', headerlines, 0);

    % trying fscanf, but it doesn't read correctly (and is not much faster)
%     fID = fopen(fullfile(vic_out_dir, wbnames(k).name));
%     % skip the headerlines
%     for kk=1:headerlines
%         fgetl(fID);
%     end
%     wb_out_all(:,:,k) = fscanf(fID, formatspec, sizeA);
%     fclose(fID);

    if mod(k,1e2) == 0
        disp(k)
%         toc
    end
    
end

return