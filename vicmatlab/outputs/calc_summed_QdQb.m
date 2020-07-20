% Calculate summed runoff and baseflow
%
% wbdir = directory containing water balance output files 
function [sQ, sQd, sQb, timevector] = calc_summed_QdQb(wbdir, A_basin)

fnames = dir([wbdir '/wb_*']);
ncells = length(fnames);
dat = readmatrix(fullfile(wbdir,fnames(1).name));
timevector = datetime(dat(:,1), dat(:,2), dat(:,3));

nt = length(timevector);
Qb = zeros(nt,ncells);
Qd = zeros(nt,ncells);

for k=1:ncells
    dat = readmatrix(fullfile(wbdir, fnames(k).name));
    Qb(:,k) = dat(:,6);
    Qd(:,k) = dat(:,5);
end

% Convert from mm/day to m^3/s
ncells = 7833;
A_cell = A_basin/ncells; % km^2, average grid cell area
conversion_factor = A_cell*1000/(24*3600);
Qb = Qb.*conversion_factor; % m^3/s for each grid cell
Qd = Qd.*conversion_factor; % m^3/s for each grid cell

sQb = sum(Qb,2);
sQd = sum(Qd,2);
sQ = sQb + sQd;

% [tm, Qm] = daily_to_monthly(timevector, sQ, 'mean')

return