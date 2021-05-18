% Calculates NSE

function NSE = myNSE(y, y_hat)

% remove NaNs
A = [y, y_hat];
i1 = isnan(A);
i2 = find(sum(i1, 2) > 0);
A(i2,:) = [];
y = A(:,1);
y_hat = A(:,2);

SSE = (y_hat - y)'*(y_hat - y);
y_mean = mean(y);
sum2 = (y - y_mean)'*(y- y_mean);
NSE = 1 - SSE/sum2;

return
