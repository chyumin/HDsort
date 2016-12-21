function [m ia ib ic] = max3D(A)
% This function expands the functionality of max.
% It returns the maximal element of A as well as its row and column index
% (first dimension and second dimension in this order).
% [m ia ib] = max2D(A)

[sa, sb, sc] = size(A);

[m0, ia, ib] = hdsort.util.max2D(A(:,:,1));
m = m0; ic = 1;
for ii = 2:sc
    [m] = hdsort.util.max2D(A(:,:,ii));
    
    if m > m0
        m0 = m;
        ic = ii;
    end
end
[m, ia, ib] = hdsort.util.max2D(A(:,:,ic));

end