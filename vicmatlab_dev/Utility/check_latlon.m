% Makes sure that lat is in decreasing order, and lon is in increasing
% order
% Not super-robust, but should work for my purposes

function result = check_latlon(lat, lon)

if lat(end) < lat(1) && lon(end) > lon(1)
    result = 0; % pass
else
    warning('Lat and lon need to be in a different order')
    result = 1; % fail
end

return