% Calculates RMSE

function RMSE = myRMSE(y, y_hat)

% remove NaNs
A = [y, y_hat];
A(isnan(y),:) = [];
A(isnan(y_hat),:) = [];
y = A(:,1);
y_hat = A(:,2);

SSE = (y_hat - y)'*(y_hat - y);
N = length(y);
RMSE = sqrt(SSE/N);

return
