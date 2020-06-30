% Aggregate forcing data
%
% Aggregates VIC meteorological forcing data from finer to coarser temporal
% resolution
%
% Example
% forcdir = '/media/jschap/HD_ExFAT/ucrb/disagg_forcings_L13';
% delt_in = 1; % hours
% delt_out = 24; % hours
% starttime = datetime(1980, 1, 1, 0, 0, 0);
% outdir = '/media/jschap/HD_ExFAT/ucrb/daily_forcings_L13';
% precip_col = 1; % need to specify because precip is summed, not averaged,
% to keep units consistant (mm/timestep)
% hpc = 1; % to use parallel processing

function [aggf, aggt] = aggregate_forcing(forcdir, delt_in, delt_out, starttime, outdir, precip_col, hpc)

% Checks
if mod(delt_out, delt_in) ~= 0 
    error('delt_out must be a multiple of delt_in')
end

fnames = dir([forcdir '/full_data_*']);
ncells = length(fnames);
f = load(fullfile(forcdir, fnames(1).name));

[nt_in, nvars] = size(f);
nt_out = nt_in*delt_in/delt_out;

% Returns error for now. Eventually should incorporate chopping as nec.
if mod(nt_out, 1) ~= 0
    error('Non-integer number of output time steps. Need to chop or pad.')
    % Find the largest number<= nt_out that is evenly divisible by
    % delt_out/delt_in
%     floor(nt_out*delt_in/delt_out);
end

endtime = starttime + hours(delt_out)*nt_out - hours(delt_out);
aggt = starttime:hours(delt_out):endtime;

mkdir(outdir)
disp(['Created directory to store outputs: ' outdir])

if hpc == 1
    c = parcluster('local');
    nw = c.NumWorkers;
    parpool(nw-1)
    parfor k=1:ncells
        f = load(fullfile(forcdir, fnames(k).name));
        f_reshape = reshape(f, delt_out/delt_in, nt_out, nvars);
        aggf = squeeze(mean(f_reshape, 1));
        aggf(:, precip_col) = sum(f_reshape(:, :, precip_col), 1);
        plotflag = 0;
        if plotflag && k == 1
            figure
            plot(aggt, aggf(:, precip_col))
            xlabel('time')
            ylabel('precip')
        end
        dlmwrite(fullfile(outdir, fnames(k).name), aggf, ' ');
    end
    aggf = zeros(length(aggt),1);
else
    for k=1:ncells
        f = load(fullfile(forcdir, fnames(k).name));
        f_reshape = reshape(f, delt_out/delt_in, nt_out, nvars);
        aggf = squeeze(mean(f_reshape, 1));
        aggf(:, precip_col) = sum(f_reshape(:, :, precip_col), 1);
        plotflag = 0;
        if plotflag && k == 1
            figure
            plot(aggt, aggf(:, precip_col))
            xlabel('time')
            ylabel('precip')
        end
        dlmwrite(fullfile(outdir, fnames(k).name), aggf, ' ');
    end
end

return