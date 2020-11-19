% Calculate summed runoff and baseflow
%
% wbdir = directory containing water balance output files 
% I've checked this code pretty thoroughly. I'm convinced it gives answers
% in the correct units - 8/10/2020 JRS
%
% Validated with the Yampa basin model - 8/12/2020 JRS

function [sQ, sQd, sQb, timevector] = calc_summed_QdQb(wbdir, prefix, A_basin, varargin)

numvarargs = length(varargin);
if numvarargs > 2
    error('The max number of optional arguments is 5')
end

optargs = {6, 7};
optargs(1:numvarargs) = varargin;
[runoff_col, baseflow_col] = optargs{:};

if runoff_col == 6 && baseflow_col == 7
    disp('Assuming runoff and baseflow are in columns 6 and 7 of outputs')
end

fnames = dir([wbdir '/' prefix '*']);

if strcmp(fnames(1).name(end-3:end), '.txt')
    ncells = length(fnames);
    dat = readmatrix(fullfile(wbdir,fnames(1).name));
    timevector = datetime(dat(:,1), dat(:,2), dat(:,3));
    nt = length(timevector);
    Qb = zeros(nt,ncells);
    Qd = zeros(nt,ncells);    
    for k=1:ncells
        dat = readmatrix(fullfile(wbdir, fnames(k).name));
        Qb(:,k) = dat(:,baseflow_col); % baseflow
        Qd(:,k) = dat(:,runoff_col); % runoff
    end
else
    ncells = length(fnames);
    dat = dlmread(fullfile(wbdir,fnames(1).name));
    timevector = datetime(dat(:,1), dat(:,2), dat(:,3));
    nt = length(timevector);
    Qb = zeros(nt,ncells);
    Qd = zeros(nt,ncells);    
    for k=1:ncells
        dat = dlmread(fullfile(wbdir, fnames(k).name));
        Qb(:,k) = dat(:,baseflow_col); % baseflow
        Qd(:,k) = dat(:,runoff_col); % runoff
    end
end

% Convert from mm/day to m^3/s
% ncells = 7833;
A_cell = A_basin/ncells; % km^2, average grid cell area
conversion_factor = A_cell*1000/(24*3600);
Qb = Qb.*conversion_factor; % m^3/s for each grid cell
Qd = Qd.*conversion_factor; % m^3/s for each grid cell

sQb = sum(Qb,2);
sQd = sum(Qd,2);
sQ = sQb + sQd;

% [tm, Qm] = daily_to_monthly(timevector, sQ, 'mean')

return