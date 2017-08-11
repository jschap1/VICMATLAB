% Plots map of basin, flow network, flow direction file, and gauge locations

% Written by Dongyue Li
% Modified by Jacob Schaperow, Aug. 11, 2017

clear, clc
cd /Users/jschapMac/Desktop/Tuolumne2/RoutInputs/

%% INPUTS

fdir=dlmread('fdir_nohead.txt'); % flow direction file, with the header removed

gage = [-121.155, 37.595]; % gage locations

% Copy these indices from the stnloc file
ind = [14, 10]; % row, col of the flow direction file where the gages are located

res=1/16; % VIC modeling resolution

% Get lat/lon of basin mask (only the pixels whose values are 1)
[mask,R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne/RoutingInputs/basinmask_coarse.asc');
[nrow, ncol] = size(mask);
[minlat, minlon] = pix2latlon(R,nrow,1);
[maxlat, maxlon] = pix2latlon(R,1,ncol);

lat = maxlat:-res:minlat;
% lat = minlat:res:maxlat;
lon = minlon:res:maxlon;

[Domain.X,Domain.Y] = meshgrid(lon,lat);

% Basin polygon outline
tuolumne = shaperead('/Users/jschapMac/Desktop/Tuolumne/Shapefiles/upper_tuolumne_wgs.shp'); 
% note that the basin boundary shapefile does not match the VIC domain 
% perfectly because of the coarse resolution

%% Plotting

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
  
% If this block is plots the basin upside-down, try reversing the order of
% the Domain.Y vector.
interval=res/2;
figure
for i=1:size(Domain.X,1)
    for j=1:size(Domain.X,2)
        % draws grid
        if fdir(i,j)>=0 % condition for being in the basin boundary
            plot([Domain.X(i,j)-interval,Domain.X(i,j)+interval,...
                Domain.X(i,j)+interval,Domain.X(i,j)-interval,Domain.X(i,j)-interval],...
                [Domain.Y(i,j)-interval,Domain.Y(i,j)-interval,...
                Domain.Y(i,j)+interval,Domain.Y(i,j)+interval,Domain.Y(i,j)-interval],'k-')            
            hold on
 
        end
           
    end
end

plot(tuolumne.X,tuolumne.Y,'k-')

quiver(Domain.X,Domain.Y,u,v)
hold on
for i=1:size(ind,1) % Plot gage locations
    % The red dot is not where I would expect it
    plot(Domain.X(1,ind(i,1)),Domain.Y(nrow-ind(i,2)+1),'r*')
    hold on
    plot(gage(i,1),gage(i,2),'bo','linewidth',2)
    hold on
%     pause
end
hold off