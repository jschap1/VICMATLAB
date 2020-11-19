% mm/day to m^3/s
%
% Converts VIC outputs (e.g. runoff and baseflow) to volumetric flowrate
%
% INPUTS
% Q = runoff depth (basin or pixel), in mm/day
% A = grid cell area, in km^2
%
% OUTPUTS
% Qv = volumetric flowrate

function Qv = mmday2cms(Q, A_cell)

cf = 1000*A_cell/(24*3600);
Qv = cf*Q;

return