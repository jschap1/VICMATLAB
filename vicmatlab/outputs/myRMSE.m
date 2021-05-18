% Calculates RMSE

function RMSE = myRMSE(y, y_hat)

% remove NaNs
A = [y, y_hat];
i1 = isnan(A);
i2 = find(sum(i1, 2) > 0);
A(i2,:) = [];
y = A(:,1);
y_hat = A(:,2);

SSE = (y_hat - y)'*(y_hat - y);
N = length(y);
RMSE = sqrt(SSE/N);

return
