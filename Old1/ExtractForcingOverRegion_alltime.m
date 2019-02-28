function AA = ExtractForcingOverRegion_alltime()

% Finds the precise forcing data that applies to the region of interest,
% from a larger forcing dataset.
%
% INPUTS
%
%
% OUTPUTS
%
%

cd '/Users/jschapMac/Desktop/UpperKaweah';
kaweah = shaperead('upper_kaweah.shp');

% ncinfo('MetNC/prec.2007.nc')
lat = ncread('MetNC/prec.2007.nc', 'lat');
lon = ncread('MetNC/prec.2007.nc', 'lon');
time = ncread('MetNC/prec.2007.nc', 'time');
prec = ncread('MetNC/prec.2007.nc', 'prec');

% Convert lon coords if they use E/W, instead of absolute value system
if max(lon)>180
    lon = lon - 360;
elseif min(lon)<-180
    % lon = lon some other conversion
end

plotflag=1;
if plotflag
    figure, hold on
    plot(kaweah.X, kaweah.Y)
    [X,Y] = meshgrid(lon, lat);
    plot(X,Y); hold off
end

% Find the cells that are in the Kaweah basin, and extract them.
figure, hold on 
plot(kaweah.X,kaweah.Y,'LineWidth',2)

plot(X,Y, '*');
hold off

mask = NaN(922, 444);
for i=1:922
    for j=1:444
        mask(i,j) = inpolygon(lon(i),lat(j),kaweah.X,kaweah.Y);
    end
end

[ind1, ind2] = find(mask);

figure, hold on 
plot(kaweah.X,kaweah.Y,'LineWidth',2)
plot(lon(ind1), lat(ind2), '*')
hold off

ncells = length(ind1);
for i=1:ncells
    data = prec(ind1(i), ind2(i), :);
    filename = ['data_' num2str(lat(ind2(i))) '_' num2str(lon(ind1(i)))];
    dlmwrite(fullfile('ClippedForcings', [filename '.txt']), data)
end


end

%%
% 
% prec_ext = prec(ind1,ind2,:); % extracted precip for the chosen points
% 
% % Next, put this into array for each grid cell
% ncells = length(ind1);
% numforcings = 4;
% for k=1:ncells
%     filename = ['data_' num2str(lat(ind2(k))) '_' num2str(lon(ind1(k)))];
%     prec_ext(lon(ind1),lat(ind2), k)
% 
% 
% prec_record = cell(ncells,2+numforcings);
% prec_record(:,1) = lon(ind1);
% prec_record{:,2} = lat(ind2);
% prec_record{:,3} = prec_ext(1, 1, :);
% 
% data_48.1875_-120.6875
% % ntimesteps by 4
% 
% end
% 
