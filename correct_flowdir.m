% Plots map of basin, flow network, flow direction file, and gauge locations

% Written by Dongyue Li
% Modified by Jacob Schaperow, Aug. 11, 2017

clear, clc
cd /Users/jschapMac/Desktop/Tuolumne2/RoutInputs/

%% INPUTS

fdir=dlmread('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/fdir_in_nohead.asc'); % flow direction file, with the header removed

gage = [-121.155, 37.595]; % gage locations

% Copy these indices from the stnloc file
ind = [3, 3]; % row, col of the flow direction file where the gages are located

res=0.007996332831; % modeling resolution

% Get lat/lon of basin mask (only the pixels whose values are 1)
[mask,R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/bmask_1k_wgs.asc');
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

% River network shapefile (optional) from an established database
rivs = shaperead('/Users/jschapMac/Desktop/Tuolumne5/LFP/rivs_clip_wgs/rivs_clip_wgs.shp');

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

% Plot rivers from database
for i=1:length(rivs)
    plot(rivs(i).X,rivs(i).Y,'r-')
end

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

%% Do the correction

fID1 = fopen('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/xcoords2change', 'a');
fID2 = fopen('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/ycoords2change', 'a');

[xx,yy] = ginput();

fprintf(fID1,'%7.4f \n',xx);
fprintf(fID2,'%6.4f \n',yy);

fclose(fID1);
fclose(fID2);

%%
% Load coordinates

xx = load('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/xcoords2change.txt');
yy = load('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/ycoords2change.txt');


%%
% (based on find_stnloc script)

[fd, R] = arcgridread('/Users/jschapMac/Desktop/Tuolumne5/LFP/FlowDir/fdir_in.asc');
coords = [xx,yy];
min_diff = 0.008; % resolution

% check each cell in the fdir file and see if it is closest to the stations
[nrows, ncols] = size(fd);

row = NaN(length(xx),1);
col = NaN(length(xx),1);
for pp=1:length(xx) % should vectorize this for speed
    
    gage = [xx(pp) yy(pp)];

    breakflag = 0; % needed bc it is a pain for Matlab to break out of nested loops
    for i=1:nrows
        for j=1:ncols
            [lat, lon] = pix2latlon(R,i,j);
            xdif = abs(lon - gage(:,1));
            ydif = abs(lat - gage(:,2));

            if xdif < min_diff && ydif < min_diff
                disp(['found point ' num2str(pp)])
                row(pp) = i;
                col(pp) = j;
                breakflag = 1;
                break;
            end

        end

        if breakflag
            break;
        end

    end

end

%% Change the values

fdirnew = fdir;

replacenum = 7*ones(length(xx),1);
replacenum(2) = 1;

for pp=1:length(xx)
    fdirnew(row(pp),col(pp)) = replacenum(pp);
end

%% Check that the changes are good

%% Save the modified flow direction file
