function [m ia ib] = min2D(A)
% This function expands the functionality of min.
% It returns the minimal element of A as well as its row and column index
% (first dimension and second dimension in this order).
[m_ i_] = min(A);
[m ib] = min(m_);
ia = i_(ib);
end