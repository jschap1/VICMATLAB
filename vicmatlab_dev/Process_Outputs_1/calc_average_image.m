function [avg_ts, avg_map] = calc_average_image(fluxfile, info_image)


for i = 1:length(info_image.vars)
    
    var = ncread(fluxfile, info_image.vars{i});
    
    % TODO: amend to handle parameters with a dimension for each soil layer
    if ndims(var) ~=3
        continue
    end
    
    var = permute(var, [2,1,3]);
    
    avg_map.(info_image.vars{i}) = mean(var, 3);
    avg_ts.(info_image.vars{i}) = squeeze(mean(mean(var,1),2));
    
end

return