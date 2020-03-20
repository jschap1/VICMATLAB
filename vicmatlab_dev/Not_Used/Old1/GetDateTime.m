function x = GetDateTime(fluxresults)

% Extracts the date time vector from VIC results.
%
% INPUTS
% VIC flux results, as from LoadVICResults. The first three columns must be
% year, month, day.
%
% OUTPUTS
% x = datetime vector for the VIC flux results

Y = fluxresults(:,1,1);
M = fluxresults(:,2,1);
D = fluxresults(:,3,1);
x = datetime(Y,M,D);

end