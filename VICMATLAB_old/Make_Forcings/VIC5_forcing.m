%% Create forcing input file

clc
clear

% ##########################################################################
file=dir('/Users/dongyueli/Desktop/fd/diss/full*'); % read in hourly disaggregated forcing
for i=1:length(file)
    lat_rec(i,1)=str2num(file(i).name(11:18));    
    lon_rec(i,1)=str2num(file(i).name(20:29));   
end

min_lat=min(lat_rec);
max_lat=max(lat_rec);
min_lon=min(lon_rec);
max_lon=max(lon_rec);

lat=([min_lat:0.0625:max_lat])';
lon=([min_lon:0.0625:max_lon])';

% get the lon/lat index of each cell in the soil file
for i=1:length(file)
    lon_ind(i,1)=find(lon==lon_rec(i,1));
    lat_ind(i,1)=find(lat==lat_rec(i,1));
end

% sp_data=dlmread(['./disaggregated_forcing/',file(1).name]);
sp_data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(1).name]);
time=[0:1:(size(sp_data,1)-1)]';

mask=nan(length(lon),length(lat));
for i=1:length(file)
    mask(lon_ind(i,1),lat_ind(i,1))=1; 
end
% ##########################################################################


for jj=1:1
    
nccreate('Sierra_Livneh_hr_19792011.nc','time',...
    'Datatype','int32',...
    'Dimensions',{'time',length(time)},...
          'Format','netcdf4_classic')

nccreate('Sierra_Livneh_hr_19792011.nc','lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','netcdf4_classic')
      
nccreate('Sierra_Livneh_hr_19792011.nc','lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','netcdf4_classic')    
     
nccreate('Sierra_Livneh_hr_19792011.nc','mask',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      
nccreate('Sierra_Livneh_hr_19792011.nc','prcp',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')      
      
nccreate('Sierra_Livneh_hr_19792011.nc','tas',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')       
      
nccreate('Sierra_Livneh_hr_19792011.nc','dswrf',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   

nccreate('Sierra_Livneh_hr_19792011.nc','dlwrf',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   

nccreate('Sierra_Livneh_hr_19792011.nc','pres',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   

nccreate('Sierra_Livneh_hr_19792011.nc','vp',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   

nccreate('Sierra_Livneh_hr_19792011.nc','wind',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
     
end

ncwrite('Sierra_Livneh_hr_19792011.nc','time',time);
ncwrite('Sierra_Livneh_hr_19792011.nc','lon',lon);
ncwrite('Sierra_Livneh_hr_19792011.nc','lat',lat);
ncwrite('Sierra_Livneh_hr_19792011.nc','mask',mask);


prcp=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   prcp(lon_ind(i,1),lat_ind(i,1),:)=data(:,1);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','prcp',prcp);
clear prcp

tas=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   tas(lon_ind(i,1),lat_ind(i,1),:)=data(:,2);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','tas',tas);
clear tas

dswrf=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   dswrf(lon_ind(i,1),lat_ind(i,1),:)=data(:,3);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','dswrf',dswrf);
clear dswrf

dlwrf=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   dlwrf(lon_ind(i,1),lat_ind(i,1),:)=data(:,4);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','dlwrf',dlwrf);
clear dlwrf

pres=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   pres(lon_ind(i,1),lat_ind(i,1),:)=data(:,6);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','pres',pres);
clear pres

vp=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   vp(lon_ind(i,1),lat_ind(i,1),:)=data(:,7);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','vp',vp);
clear vp


wind=nan(length(lon),length(lat),length(time));
tic
for i=1:length(file)
   data=[];
   data=dlmread(['/Users/dongyueli/Desktop/fd/diss/',file(i).name]);
   wind(lon_ind(i,1),lat_ind(i,1),:)=data(:,8);   
end
toc
ncwrite('Sierra_Livneh_hr_19792011.nc','wind',wind);
clear wind



for jj=1:1 % generate a bulk nc file that has hourly data from all the years

    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'time','units','hours since 1979-01-01');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'time','calendar','proleptic_gregorian');


    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lon','standard_name','longitude');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lon','long_name','longitude of grid cell center');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lon','units','degrees_east');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lon','axis','X');

    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lat','standard_name','latitude');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lat','long_name','latitude of grid cell center');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lat','units','degrees_north');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'lat','axis','Y');

    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'mask','_FillValue',9.969209968386869e+36);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'mask','comment','NaN indicates cell is not active');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'mask','long_name','fraction of grid cell that is activedomain mask');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'mask','note','unitlessunitless');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'prcp','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'prcp','long_name','PREC');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'prcp','column',0);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'prcp','units','mm/step');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'prcp','description','PREC');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'tas','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'tas','long_name','AIR_TEMP');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'tas','column',1);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'tas','units','C');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'tas','description','AIR_TEMP');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'dswrf','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dswrf','long_name','SWDOWN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dswrf','column',2);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dswrf','units','S_/m2');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dswrf','description','SWDOWN');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'dlwrf','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dlwrf','long_name','LWDOWN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dlwrf','column',3);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dlwrf','units','S_/m2');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'dlwrf','description','LWDOWN');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'pres','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'pres','long_name','PRESSURE');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'pres','column',5);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'pres','units','KPa');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'pres','description','PRESSURE');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'vp','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'vp','long_name','VP');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'vp','column',6);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'vp','units','KPa');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'vp','description','VP');

    % ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
    %     'wind','_FillValue','NaN');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'wind','long_name','WIND');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'wind','column',7);
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'wind','units','m/s');
    ncwriteatt('Sierra_Livneh_hr_19792011.nc',...
        'wind','description','WIND');

end



%% Read in all the forcing and disaggregate to each year as VIC5 forcing

clc
clear


doyr=([365 366 365 365 365 366 365 365 365 366 365 365 365 366 365 365 365 ...
     366 365 365 365 366 365 365 365 366 365 365 365 366 365 365 365])';
yr=[1979:1:2011];

% time=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc','time');
lon=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc','lon');
lat=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc','lat');
mask=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_1979.nc','mask');



for i=1:length(yr)

    tic
    
    year=num2str(yr(i));
    
    time=[0:3:24*doyr(i)-2]'; % generate 3hr forcing for each year
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','prcp');
    prcp=data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0);
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','tas');
    tas=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','dswrf');
    dswrf=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','dlwrf');
    dlwrf=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','pres');
    pres=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','vp');
    vp=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    
    data=ncread('/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_AllYears.nc','wind');
    wind=(data(:,:,24*sum(doyr(1:i-1))+1:3:24*sum(doyr(1:i))-2)+data(:,:,24*sum(doyr(1:i-1))+2:3:24*sum(doyr(1:i))-1)+...
        data(:,:,24*sum(doyr(1:i-1))+3:3:24*sum(doyr(1:i))-0))/3;
    clear data
    

    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'time',...
    'Datatype','int32',...
    'Dimensions',{'time',length(time)},...
          'Format','netcdf4_classic')
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'time',time);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'lon',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon)},...
          'Format','netcdf4_classic')
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'lon',lon);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'lat',...
    'Datatype','double',...
    'Dimensions',{'lat',length(lat)},...
          'Format','netcdf4_classic')
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'lat',lat);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'mask',...
    'Datatype','double',...
    'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'mask',mask);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'prcp',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')      
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'prcp',prcp);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'tas',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')       
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'tas',tas);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'dswrf',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'dswrf',dswrf);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'dlwrf',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'dlwrf',dlwrf);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'pres',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'pres',pres);
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'vp',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'vp',vp)
      
    nccreate(['Sierra_Livneh_hr_2011_',year,'.nc'],'wind',...
    'Datatype','single',...
    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',length(time)},...
          'Format','netcdf4_classic')   
      ncwrite(['Sierra_Livneh_hr_2011_',year,'.nc'],'wind',wind)
    
    clear prcp tas dswrf dlwrf pres vp wind time
     
    
    
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'time','units',['hours since ',year,'-01-01']); % the date here must be right
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'time','calendar','proleptic_gregorian');


    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lon','standard_name','longitude');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lon','long_name','longitude of grid cell center');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lon','units','degrees_east');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lon','axis','X');

    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lat','standard_name','latitude');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lat','long_name','latitude of grid cell center');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lat','units','degrees_north');
    ncwriteatt(['Sierra_Livneh_hr_2011_',year,'.nc'],...
        'lat','axis','Y');


    
    toc
    
  
end


yr=[1979:1:2011];
for i=1:length(yr)
    
year=num2str(yr(i));

ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'time','units',['hours since ',year,'-01-01']); % the date here must be right
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'time','calendar','proleptic_gregorian');


ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lon','standard_name','longitude');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lon','long_name','longitude of grid cell center');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lon','units','degrees_east');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lon','axis','X');

ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lat','standard_name','latitude');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lat','long_name','latitude of grid cell center');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lat','units','degrees_north');
ncwriteatt(['/Users/dongyueli/Desktop/VIC/setup/Livneh_sierra_VIC5/forcings/Sierra_Livneh_hr_2011_',year,'.nc'],...
    'lat','axis','Y');

end
