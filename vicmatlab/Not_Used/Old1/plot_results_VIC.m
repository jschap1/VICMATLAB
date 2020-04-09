%% 
% Load
loadfold = '/Users/jschapMac/Desktop/vis_Stehekin';
cd(loadfold)
load fluxresults % dim = timesteps x numvars x gridcells
load snowresults % dim = timesteps x numvars x gridcells
load routresults % dim = timesteps x numvars

%%


%% 
% Separate by variable
tsindex = 1; % 1 for yearly, 2 for monthly, 3 for daily
numlayers = 3; % number of soil layers

YEAR = fluxresults(:,1,:);
MONTH = fluxresults(:,2,:); 
DAY = fluxresults(:,3,:);     
OUT_PREC = fluxresults(:,4,:);        
OUT_EVAP = fluxresults(:,5,:);       
OUT_RUNOFF = fluxresults(:,6,:);  
OUT_BASEFLOW = fluxresults(:,7,:);  
OUT_WDEW = fluxresults(:,8,:);       
OUT_SOIL_LIQ = fluxresults(:,1,:); % layers 1-3 are grouped
OUT_NET_SHORT = fluxresults(:,1,:);   
OUT_R_NET = fluxresults(:,1,:);       
OUT_EVAP_CANOP = fluxresults(:,1,:);  
OUT_TRANSP_VEG = fluxresults(:,1,:);  
OUT_EVAP_BARE = fluxresults(:,1,:);   
OUT_SUB_CANOP = fluxresults(:,1,:);   
OUT_SUB_SNOW = fluxresults(:,1,:);    
OUT_AERO_RESIST = fluxresults(:,1,:);         
OUT_SURF_TEMP = fluxresults(:,1,:);   
OUT_ALBEDO = fluxresults(:,1,:);      
OUT_REL_HUMID = fluxresults(:,1,:);   
OUT_IN_LONG = fluxresults(:,1,:);     
OUT_AIR_TEMP = fluxresults(:,1,:);    
OUT_WIND = fluxresults(:,1,:);

%%
% Radiation

fluxresults


