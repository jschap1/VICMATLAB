% Array to map (VIC)
%
% INPUTS:
% Qd = [npix, nt] array of VIC outputs for each grid cell
% nanmask = basin mask (ones and NaNs)
%
% Be sure that sum(nanmask(:)==1) is equal to npix.

function Qdmap4 = array2map_vic(Qd, nanmask)

[npix, nt] = size(Qd);
[nlon, nlat] = size(nanmask);

Qdmap = zeros(nlon, nlat, nt);
for t=1:nt
    Qdmap(:,:,t) = nanmask;
end

ind = nanmask==1;

Qdmap2 = struct();
for t=1:nt
    Qdmap2.(['t' num2str(t)]) = Qdmap(:,:,t);
    Qdmap2.(['t' num2str(t)])(ind) = Qd(:,t);
end

% figure, plotraster(masklon, masklat, Qdmap2.t1, '', '', '')
    
Qdmap3 = struct2cell(Qdmap2);

Qdmap4 = zeros(nlon, nlat, nt);
for t=1:nt
    Qdmap4(:,:,t) = Qdmap3{t};
end

return