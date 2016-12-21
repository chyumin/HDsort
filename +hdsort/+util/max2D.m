function [m ia ib] = max2D(A, N)
% This function expands the functionality of max.
% It returns the maximal element of A as well as its row and column index
% (first dimension and second dimension in this order).
% [m ia ib] = max2D(A, N)

if nargin == 1
    [m_ i_] = max(A);
    [m ib] = max(m_);
    ia = i_(ib);
elseif isnumeric(N)
    [sA idxsA] = sort(A(:), 'descend');
    m = sA(1:N);
    [ia ib ] = ind2sub(size(A), idxsA(1:N));
elseif strcmp(N, 'each column')
    [m_ ia_] = max(A);
    ib_ = 1:size(A,2);
    
    [m idx] = sort(m_, 'descend');
    ia = ia_(idx);
    ib = ib_(idx);
elseif strcmp(N, 'each row')
    [m_ ib_] = max(A');
    ia_ = 1:size(A,1);
    
    [m idx] = sort(m_, 'descend');
    ia = ia_(idx);
    ib = ib_(idx);
else
    error('Wrong N!')
end

end