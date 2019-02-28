%% Plot map of basin, flow network, flow direction file, and gauge locations

cd /Users/jschapMac/Desktop/Sierra

%% INPUTS

fdir=dlmread('sierra_fcst.txt'); % flow direction file

load('sierra_flow_domain.mat'); %  
load('sierra_normal_shape.mat'); % 

gage = []; % gage locations

ind = []; % row, col of the flow direction file where the gages are located

res=1/16; % VIC modeling resolution

Domain.X
Domain.Y % meshgrid results of lat and lon vectors spanning the study area

Domain.ClipFlow

axesm sinusoid;

Domain.ClipFlow % flow direction file, with cells not in the basin set to -1

% Basin polygon outline
% (sierra.X,sierra.Y)

%% CORRECTING FLOW FILE, FIND RIGHT POSITION FOR EACH FORECAST PT
% Written by Dongyue Li, 2017

clc
clear
% copying Li.flow from txt to a=[];
% a=dlmread('/Users/dongyueli/Desktop/VIC/data_prepare/Routing/Li_flow.txt');


% Inputs are flow direction and basin mask
% load('Kostas_geo.mat');
load('sierra_flow_domain.mat'); 
load('sierra_normal_shape.mat');
 
u=zeros(size(fdir,1),size(fdir,2));
v=u;
for i=1:size(u,1)
    for j=1:size(u,2)
        
        if fdir(i,j)==1
            u(i,j)=0;
            v(i,j)=1;
        end
        
        if fdir(i,j)==2
            u(i,j)=1;
            v(i,j)=tand(45);
        end
        
        if fdir(i,j)==3
            u(i,j)=1;
            v(i,j)=0;
        end
        
        if fdir(i,j)==4
            u(i,j)=1;
            v(i,j)=-tand(45);
        end
        
        if fdir(i,j)==5
            u(i,j)=0;
            v(i,j)=-1;
        end
        
        if fdir(i,j)==6
            u(i,j)=-1;
            v(i,j)=-tand(45);
        end
        
        if fdir(i,j)==7
            u(i,j)=-1;
            v(i,j)=0;
        end
        
        if fdir(i,j)==8
            u(i,j)=-1;
            v(i,j)=tand(45);
        end
        
    end
end
 
 
gage=[
-122.186    40.289
-121.547    39.522
-121.274124 39.235172
-121.183    38.683
-121.044    38.5
-120.719    38.313
-120.637    37.852
-120.441    37.666
-120.331    37.522
-119.724312 36.984394
-119.335    36.831
-119.003    36.412
-118.922    36.061
-118.484    35.639];
 
ind=[15 80;
    27 68;
    28 62;
    32 54;
    34 51;
    40 48;
    42 42;
    44 38;
    47 36;
    56 27;
    62 24;
    67 17;
    68 12;
    75 5];
 
interval=res/2;
figure
for i=1:size(Domain.X,1)
    for j=1:size(Domain.X,2)
        
        if Domain.ClipFlow(i,j)~=-1
            plot([Domain.X(i,j)-interval,Domain.X(i,j)+interval,...
                Domain.X(i,j)+interval,Domain.X(i,j)-interval,Domain.X(i,j)-interval],...
                [Domain.Y(i,j)-interval,Domain.Y(i,j)-interval,...
                Domain.Y(i,j)+interval,Domain.Y(i,j)+interval,Domain.Y(i,j)-interval],'k-')            
            hold on
 
        end
           
    end
end
for i=1:25
    plot(sierra(i).X,sierra(i).Y,'k-')
    hold on
end
quiver(Domain.X,Domain.Y,u,v)
hold on
for i=1:size(ind,1)
    plot(Domain.X(1,ind(i,1)),Domain.Y(105-ind(i,2)+1,1),'r*')
    hold on
    plot(gage(i,1),gage(i,2),'bo','linewidth',2)
    hold on
%     pause
end
hold off