function [] = LoadVICResults()

% Load results from VIC simulation and the routing model

cd '/Users/jschapMac/Desktop/results_VIC/longrun'
fluxnames = dir('fluxes*');
ncells = length(fluxnames);
tmp = dlmread(fluxnames(1).name);  
fluxresults = NaN([size(tmp),ncells]);

for i=1:ncells
    fluxresults(:,:,i) = dlmread(fluxnames(i).name);  
end

snownames = dir('snow*');
ncells = length(snownames);
tmp = dlmread(snownames(1).name);  
snowresults = NaN([size(tmp),ncells]);

for i=1:ncells
    snowresults(:,:,i) = dlmread(snownames(i).name);  
end

%%
% Load results from routing model

cd '/Users/jschapMac/Desktop/results_rout/lrun'
prefix = 'STEHE';
units = 'mm'; % mm, cfs
timestep = 'month'; % day, month, year

switch units
    case 'cfs'
       routname = [prefix '.' timestep];
    case 'mm'
       routname = [prefix '.' timestep '_mm'];
end

routresults = dlmread(routname);

%% 
% Save
savefold = '/Users/jschapMac/Desktop/vis_Stehekin/longrun';
cd(savefold)
save fluxresults
save snowresults
save routresults

end
