% Creates input files for the Lohmann routing model
%
% Written 10/23/2019 JRS

% wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Raw_EB_FS_1980-2011/wb';
% wb_out_dir = '/Volumes/HD4/SWOTDA/Data/IRB/EB_1980-2019_CORR/out/wb';
% wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/Classic_VICGlobal/Raw/wb';
wb_out_dir = '/home/jschap/Documents/ESSD/distrib_cal/WY2002/wb';

fnames = dir([wb_out_dir '/wb_*']);

ncells = length(fnames);
% outdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/Classic_VICGlobal/Raw/rout_in';
outdir = '/home/jschap/Documents/ESSD/distrib_cal/rout_out';
mkdir(outdir)
for k=1:ncells
    
    dat = dlmread(fullfile(wb_out_dir, fnames(k).name), '\t', 3, 0);
    nt = size(dat,1);
%     A = [dat(:,1:3), zeros(nt,1), zeros(nt,1), dat(:,6), dat(:,7)];
    
    A = [dat(:,1:3), zeros(nt,1), zeros(nt,1), dat(:,5), dat(:,6)];
    outname = ['fluxes_' fnames(k).name(4:end)]; % if name is wb_
%     outname = ['fluxes_' fnames(k).name(10:end)]; % if name is fluxes_ or wb_daily
    dlmwrite(fullfile(outdir, outname), A, '\t')
    
end


(1-100/170)/(-1.5)