% cfs to cms
%
% INPUTS
% Qcfs = volumetric flowrate, English units
%
% OUTPUTS
% Qcms = volumetric flowrate, metric units

function Qcms = cfs2cms(Qcfs)

Qcms = Qcfs*(12/39.37)^3;

return