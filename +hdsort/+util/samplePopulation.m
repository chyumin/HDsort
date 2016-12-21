function [idx, y, ny, P] = samplePopulation(n, k, varargin)
% This function returns the indices of k random elements of each unique 
% element in input vector n.
% 
P.debug = false;
P.permute = false;
P = hdsort.util.parseInputs(P,varargin, 'error');     

logical_input = false;
if islogical(n)
    logical_input = true;
end

cellInput = iscell(n); m = {};
if cellInput
    m = n;
    n = [];
    inputtype = class(m{1});
    for jj = 1:numel(m)
        assert(isa(m{jj}, inputtype), 'Input cell must be of a consistent variable type!')
        
        m{jj} = m{jj}(:);
        n = [n jj*ones(1, numel(m{jj}))];
    end
end
n = n(:)';

unique_n = unique(n);
uN = numel(unique(n));

idx = [];
ny = zeros(1, uN);
for ii = 1:uN
    nIdx = find(n==unique_n(ii));
    
    if numel(nIdx) >= k
        ny(ii) = k;
        if ~P.permute
            idx = [idx sort(randsample(nIdx, k))];
        else
            idx = [idx randsample(nIdx, k)];
        end
    else
        idx = [idx nIdx];
        ny(ii) = numel(nIdx); 
    end 
end

y = n(idx);


if cellInput
    [unique_n, Nunique_n] = hdsort.util.unique(n)
    
    uy = unique(y);
    csum = [0; cumsum(Nunique_n)];
    z = [];
    idx2 = {};
    for jj = 1:numel(uy)
        idx2{jj} = idx(y==uy(jj))-csum(jj)
        z_ =  m{jj}(idx2{jj});
        z = [z; z_(:)];
    end
    y = z;
end

if logical_input
    y = logical(y);
end

end