% Removes empty (all NaN) rows or columns from a matrix
%
% 8/10/2020 JRS

function B = remove_empty_rowcol(A)

% empty_cols = find(nansum(A,1)==0);
% empty_rows = (nansum(A,2)==0);

non_empty_cols = find(nansum(A,1)~=0);
non_empty_rows = (nansum(A,2)~=0);

B = A(non_empty_rows, :);
B = B(:,non_empty_cols);

return