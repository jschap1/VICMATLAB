% Subset forcings (in time)
%
% start_date = initial time for input forcings
% t0 = initial time for subsetting
% tf = final time for subsetting

function subset_forcings_time(indir, outdir, start_date, t0, tf, prefix, time_step)

mkdir(outdir)
disp(['Created directory ', outdir, ' for outputs']);

%% Subset

forcnames = dir([fullfile(indir, prefix) '*']);

origfile = fullfile(indir, forcnames(1).name);
forc = load(origfile);
[nt, ~] = size(forc);
ncells = length(forcnames);
dt = hours(time_step);
date_vect = start_date:dt:(start_date + (nt-1)*dt);
date_vect = date_vect';
[~,start_ind] = ismember(t0, date_vect);
[~,end_ind] = ismember(tf, date_vect);

for k=1:ncells
    origfile = fullfile(indir, forcnames(k).name);
    forc = load(origfile);    
    newforc = forc(start_ind:end_ind,:);
    dlmwrite(fullfile(outdir, forcnames(k).name), newforc, ' ');
end

end
