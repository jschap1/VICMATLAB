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

if strcmp(class(basinmask), 'double')
    basinmask = double(basinmask);
end

disp('masking')

% Do the masking
output = load(savename);
nt = length(output.timevector);
[nlon, nlat] = size(basinmask);
output_map_masked = NaN(nlon, nlat, nt);
disp(['There are ' num2str(nt) ' time steps'])
for t=1:nt
    output_map = xyz2grid(output.info.lon, output.info.lat, output.var1(:,t));
    
    if t==1
        if size(output_map,1) ~= size(basinmask,1) || size(output_map,2) ~= size(basinmask,2)
            % Sometimes the basinmask has an extra row and column compared to
            % the 'output_map'. Remove it if necessary. This issue could
            % probably be avoided by being more careful about coordinate
            % systems
            disp('Removing extra rows and/or columns from basinmask')
            basinmask = remove_empty_rowcol(basinmask);
            [nlon, nlat] = size(basinmask);
            output_map_masked = NaN(nlon, nlat, nt);
        end
    end
        
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
    disp('saved mask for current variable')
end



return
