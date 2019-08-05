% Read forcing data
%
% INPUTS
% frequency = 'hourly', 'daily', or 'monthly'
%
% OUTPUTS
% forcings = structure array comprising maps of forcing data
%
% TODO: add capacity to compute daily and monthly averages, just reading in
% some of the data. Use fgetl creatively to do this.

function forcings = readforc(forcenames, lon, lat, resolution, frequency)

ncells = length(forcenames);
forc_now = NaN(ncells, 7);

for k=1:ncells

    fID = fopen(forcenames(k).name);
    tline = fgetl(fID);
    forc_now(k,:) = str2num(tline); % forcings for the current time step
    fclose(fID);
    
    if mod(k, 1e3) == 0
        disp(k)
    end

end
% takes 8.5 minutes for one time step, 76800 cells, uses 4 MB of RAM

temp_mat = xyz2grid(lon, lat, forc_now(:,1));
prec_mat = xyz2grid(lon, lat, forc_now(:,2));
ps_mat = xyz2grid(lon, lat, forc_now(:,3));
sw_mat = xyz2grid(lon, lat, forc_now(:,4));
lw_mat = xyz2grid(lon, lat, forc_now(:,5));
vp_mat = xyz2grid(lon, lat, forc_now(:,6));
wind_mat = xyz2grid(lon, lat, forc_now(:,7));
Rmat1 = makerefmat(min(lon), min(lat), resolution, resolution);
[tlon, tlat] = pixcenters(Rmat1, size(temp_mat));

forcings.temp = temp_mat;
forcings.prec = prec_mat;
forcings.ps = ps_mat;
forcings.sw = sw_mat;
forcings.lw = lw_mat;
forcings.vp = vp_mat;
forcings.wind = wind_mat;
forcings.R = Rmat1;
forcings.lon = tlon;
forcings.lat = tlat;

return