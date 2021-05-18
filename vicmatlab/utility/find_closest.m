% Find closest
%
% Given a coordinate (x,y), finds the closest coordinate to (x,y) in a list
% of coordinates, given a search range 
%

function [sind, mindist] = find_closest(coord, coordlist, search_range)

lon1 = coord(1);
lat1 = coord(2);

lon = coordlist(:,1);
lat = coordlist(:,2);

sind = find((abs(lon - lon1) <= search_range/2) & (abs(lat - lat1) <= search_range/2));
switch length(sind)
    case 0
        error('no nearby points found. use a larger search radius')
    case 1
        disp('one nearby point found')
        mindist = sqrt((lon(sind) - lon1).^2 + (lat(sind) - lat1).^2);
    otherwise
        disp('multiple nearby points found')
        distances = sqrt((lon(sind) - lon1).^2 + (lat(sind) - lat1).^2);
        [mindist, minind] = min(distances);
        sind = sind(minind); % keep the closest station
end

return