% Rescales Y to match the mean of X
%
% Keeps the spatial pattern of Y

function Yprime = rescale_mean(X, Y, varargin)

if nargin > 3
    error('maximum number of parameters is 3')
end

optargs = {1};
optargs(1) = varargin;
[negative_allowed] = optargs{:};

if negative_allowed
    b = mean(Y) - mean(X);
    Yprime = Y - b;
else
    t = min(Y) + 1e-3;
    s = (mean(X) - t)/mean(Y);
    Yprime = s*Y + t;
end

return
