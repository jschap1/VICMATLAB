%% Remove headers

% VIC 5 Classic output files come with three lines of header rows. In order
% to run the VIC routing model with VIC 5 outputs, these header rows need
% to be removed.
%
% Also useful:
% Remove extension from VIC output files
% for file in fluxes*.txt; do mv "$file" "${file%.txt}"; done

function outdir = remove_headers(indir, outdir, prefix)

orignames = dir(fullfile(indir, [prefix '*']));
ncells = length(orignames);
header_rows = 3;
for k=1:ncells
    data = dlmread(fullfile(indir, orignames(k).name), '\t', header_rows, 0);
    dlmwrite(fullfile(outdir, orignames(k).name), data, ' ');
    if mod(k, 1000)==0
        disp(['Processed ' num2str(k) ' of ' num2str(ncells) ' files'])
    end
end