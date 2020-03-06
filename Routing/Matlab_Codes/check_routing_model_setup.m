% May 6, 2019 JRS
% Check routing model setup before running
% 
% In development, supercedes a previous version

% flow directions file
[fd, R] = arcgridread('./Data/IRB/ROUT/irb.flowdir.asc');

% station locations file
stations_name = './Data/IRB/ROUT/irb.stnloc';
fID = fopen(stations_name, 'r');
ll = fgetl(fID);
while ~isempty(ll)
    if strcmp(ll, 'NONE')
        break;
    end
    ll = fgetl(fID);
end