% Find line
%
% July 17, 2020
% Finds the line in a text file containing the phrase (findstring)
% Note: only finds the first instance

function i = find_line(fname, findstring)

A = read_global_param_file(fname);

nrows = length(A);
for i=1:(nrows-1) % avoids issues with last line
    temp = strsplit(A{i}, findstring);
    if length(temp) > 1
        break
    end
end

return