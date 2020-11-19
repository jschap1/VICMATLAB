% ParLoad
%
% Loads VIC outputs using a parallel pool
% Enter n = 0 if you don't want to use a parallel pool

function out1 = parload(outdir, prefix, varargin)

    headerrows = 0;

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
    elseif n==0
        disp(['n = ' num2str(n)])
        disp('Not using a parallel pool')
    end
        
    wbnames = dir(fullfile(outdir, [prefix '*']));
    ncells = length(wbnames);
    out0 = dlmread(fullfile(outdir, wbnames(1).name), '\t', headerrows, 0);
    if size(out0,2) == 1
        out0 = dlmread(fullfile(outdir, wbnames(1).name), ' ', headerrows, 0);
    end
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
    
    if n>0
        
        % Check if parallel pool exists
        if isempty(gcp('nocreate'))
            parpool(n);
        end
        % delete(gcp('nocreate'))

        parfor k=1:ncells
           out1(:,k,:) = dlmread(fullfile(outdir, wbnames(k).name), '\t', [headerrows, C1, nt+(headerrows - 1), C2]);
        end
    
    else
        
        try
        for k=1:ncells
            out1(:,k,:) = dlmread(fullfile(outdir, wbnames(k).name), '\t', [headerrows, C1, nt+(headerrows - 1), C2]);
        end
        catch
        for k=1:ncells
            out1(:,k,:) = dlmread(fullfile(outdir, wbnames(k).name), ' ', [headerrows, C1, nt+(headerrows - 1), C2]);
        end            
        end
        
    end
end