function [U nU ia ib] = unique(A)
% This function expands the functionality of unique.
% It returns as a second output the number of occurences of the elements
% specified in the first output.

[U ia ib] = unique(A(:));

nU = zeros(length(U), 1);
for ii = 1:length(U)
    nU(ii) = sum(A(:)==U(ii));
end

end