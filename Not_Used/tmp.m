% For profiling 

ncells = 10;
for x=1:ncells
    savename1 = fullfile(forcdir, ['Forcings_' num2str(lonlat(x,2), precision) '_' num2str(lonlat(x,1), precision) '.txt']);
    forcing_out = zeros(ndays*24,7);
    reset(ds_temp);
    reset(ds_prec);
    reset(ds_ps);
    reset(ds_sw);
    reset(ds_lw);
    reset(ds_vp);
    reset(ds_wind);
    for d=1:ndays
        TEMP = read(ds_temp);
        PREC = read(ds_prec);
        PS = read(ds_ps);
        SW = read(ds_sw);
        LW = read(ds_lw);
        VP = read(ds_vp);
        WIND = read(ds_wind);

        forcing_out(d1(d):d2(d),1) = table2array(TEMP(x,:))'; % temperature
        forcing_out(d1(d):d2(d),2) = table2array(PREC(x,:))'; % precip
        forcing_out(d1(d):d2(d),3) = table2array(PS(x,:))'; % air pressure
        forcing_out(d1(d):d2(d),4) = table2array(SW(x,:))'; % shortwave
        forcing_out(d1(d):d2(d),5) = table2array(LW(x,:))'; % longwave
        forcing_out(d1(d):d2(d),6) = table2array(VP(x,:))'; % vapor pressure
        forcing_out(d1(d):d2(d),7) = table2array(WIND(x,:))'; % wind speed
    end
    dlmwrite(savename1, forcing_out)
end

