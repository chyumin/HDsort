
function [i, j, val] = matrixArgMax(M)
    % Computes the argmax in (row, col) of a matrix
    [m, i_] = max(M,[],1);
    [val, j] = max(m,[],2);
    val = squeeze(val);
    j   = squeeze(j);
    i = zeros(size(i_,3),1);
    for k=1:size(i_,3)
        i(k)=i_(1,j(k),k);
    end