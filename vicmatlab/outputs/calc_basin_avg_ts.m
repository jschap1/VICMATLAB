% Basin average time series
%
% An accurate way to calculate basin-average time series from a time series
% of basin maps

function ts = calc_basin_avg_ts(mapts, calcsum)

[nlon, nlat, nt] = size(mapts);

temporary = reshape(mapts, nlon*nlat, nt);
ts = nanmean(temporary, 1)'; % This method calculates the average totally accurately...

if calcsum
    ts = nansum(temporary,1)';
end

return