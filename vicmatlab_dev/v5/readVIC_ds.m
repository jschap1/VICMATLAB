% Read VIC tall
%
% Reads VIC output files, taking advantage of Matlab's big data functions
% Performs spatial and temporal averaging

function [avg_map, avg_ts, sum_ts, t, ds] = readVIC_ds(indir, nvars, ncells, nt)

ds = tabularTextDatastore(indir);

t = tall(ds);

ds.ReadSize = 'file';

% Average maps
reset(ds)
avg_map = zeros(ncells, nvars);

% Since our license was renewed, I can't call read() from within parfor. 
% It seems like the license works but not in a parfor loop.
% parfor k = 1:ncells
for k = 1:ncells
    T = read(ds);
%     if any(isnan(table2array(T)))
%         pause;
%         disp('pause')
%         1;
%     end
    avg_map(k, :) = nanmean(table2array(T));
end
% However, it does seem to work in the loop below...
%
% Handle NaNs (this handling makes speckles appear in the output by
% completely ignoring a grid cell with any NaN values - I think)
%
% It would be better to use as much information as is available to
% avoid having empty pixels.

% Average time series
reset(ds)
sum_ts = zeros(nt, nvars);
count_missing = 0;
for k = 1:ncells
% parfor k = 1:ncells
    T = read(ds);
    if ~isempty(T)
        newvals = table2array(T);
        
        if any(isnan(newvals(:)))
            count_missing = count_missing + 1;
            continue;
        end
        
        sum_ts = sum_ts + newvals;
    else
        count_missing = count_missing + 1;
    end
end
avg_ts = sum_ts/(ncells-count_missing);
disp(['There were ' num2str(count_missing) ' missing days for the temporal average']);

% save('/Volumes/HD4/SWOTDA/Data/UMRB/Classic_Livneh_met_L15/Raw/Processed/readVIC_outputs.mat', 'avg_map', 'avg_ts', 'sum_ts', 't', 'ds')

return