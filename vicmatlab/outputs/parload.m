% ParLoad
%
% Loads VIC outputs using a parallel pool
%
% Should update to only load specific columns

function out1 = parload(outdir, prefix, varargin)

    numvarargs = length(varargin);
    if numvarargs > 2
        error('The max number of optional arguments is 2')
    end

    optargs = {2, []};
    optargs(1:numvarargs) = varargin;
    [n, col] = optargs{:};
    
    disp(['Using n = ' num2str(n)])

    if n>6
        disp(['n = ' num2str(n)])
        disp('Consider using fewer processors')
    end
    
    % Check if parallel pool exists
    if isempty(gcp('nocreate'))
        parpool(n);
    end
    % delete(gcp('nocreate'))
    
    wbnames = dir(fullfile(outdir, [prefix '*']));
    ncells = length(wbnames);
    out0 = dlmread(fullfile(outdir, wbnames(1).name), '\t', 3, 0);
    [nt, nvar] = size(out0);
    
    if isempty(col)
        % Read all columns
        C1 = 0;
        C2 = nvar-1;
        out1 = zeros(nt, ncells, nvar);
    else
        % Read a single column
        C1 = col-1;
        C2 = col-1;
        out1 = zeros(nt, ncells);
    end
    
    
    parfor k=1:ncells
       out1(:,k,:) = dlmread(fullfile(outdir, wbnames(k).name), '\t', [3, C1, nt+2, C2]);
    end
    
end