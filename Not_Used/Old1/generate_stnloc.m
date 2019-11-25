%% GENERATING STATION FILE FOR GAGES
% Written by Dongyue Li, 2017

% I HAVE NOT BEEN ABLE TO DO ANYTHING USEFUL WITH THIS (JACOB)

clc
clear

% Enter coordinates of gages
gage=[-121.151    37.6];





% Load basin mask
% lon = load('lon.mat');
% lat = load('lat.mat');

for i=1:size(gage,1)
    
    col_diff=[];
    row_diff=[];
    col_diff=abs(gage(i,1)-lon);
    row_diff=abs(gage(i,2)-lat);
    
    ind(i,1)=find(col_diff==min(col_diff));
    
    a=[];
    a=find(row_diff==min(row_diff));
    
    if(length(a)~=1)
        row_diff(a(2,1))=999;
    end
    ind(i,2)=length(row_diff)-find(row_diff==min(row_diff))+1;
    
end
 
for i=1:size(ind,1)
    gage_no(i,1)=i+100;
    gage_flow(i,1)=-9000-i;    
end
    
fid=fopen('/Users/jschapMac/Desktop/fp.txt','w'); % change .txt to ascii file in the end
for i=1:size(ind,1)
    fprintf(fid,'%1i %5s             %3i %3i %5i\n',...
        1,['GS',num2str(gage_no(i))],ind(i,1),ind(i,2),gage_flow(i,1));
    fprintf(fid,'%4s\n','NONE');
end
fclose(fid);