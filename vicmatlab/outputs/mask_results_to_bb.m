% Mask results to basin boundary
%
% INPUTS
% Basin mask
% Name of .mat file where load_vic_output saved the processed VIC outputs
% Optional: filename to save the cropped VIC output 
%
%
% If the VIC results are for more grid cells than the basin covers, use
% this function to mask the VIC results to the basin boundary

function [output_map_masked, R] = mask_results_to_bb(savename, basinmask, varargin)

% Do the masking
output = load(savename);
nt = length(output.timevector);
[nlon, nlat] = size(basinmask);
output_map_masked = NaN(nlon, nlat, nt);
disp(['There are ' num2str(nt) ' time steps'])
for t=1:nt
    output_map = xyz2grid(output.info.lon, output.info.lat, output.var1(:,t));
    output_map_masked(:,:,t) = output_map.*basinmask;
    if mod(t,1000)==0
        disp(t)
    end
end

% Make georeferencing matrix
resolution = max(abs(diff(output.info.lat)));
R = makerefmat(min(output.info.lon), min(output.info.lat), resolution, resolution);

% geotiffwrite('./evap.tif', output_map_masked(:,:,1), R);

% Save the output (optional)
if ~isempty(varargin)
    save(varargin{1}, 'output_map_masked', 'R')
end

return
