% Calculates KGE

function KGE = myKGE(y, y_hat)

% remove NaNs
A = [y, y_hat];
i1 = isnan(A);
i2 = find(sum(i1, 2) > 0);
A(i2,:) = [];
y = A(:,1);
y_hat = A(:,2);

r = corr(y,y_hat); % correlation
alpha1 = std(y_hat)/std(y); % relative variability
beta1 = mean(y_hat)/mean(y); % bias
KGE = 1 - sqrt((r - 1)^2 + (alpha1 - 1)^2 + (beta1 - 1)^2);

return
