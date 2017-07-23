% Reads in Stehekin basin forcings for visualization

multipliers = [40,100,100,100];
filename = '/Users/jschapMac/Desktop/Stehekin_4_2/forcing/data_48.1875_-120.6875';

forcings = binread(filename, [41638,4]);
for i=1:size(forcings,2)
    forcings(:,i) = forcings(:,i).*multipliers(i); % should i multiply or divide by the multiplier? Divide seems more likely.
end
varnames = {'PREC','TMAX','TMIN','WIND'}';
FORCINGS = struct('forcings', forcings, 'varnames', varnames);


