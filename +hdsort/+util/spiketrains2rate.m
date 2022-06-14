function [rate, edges] = spiketrain2rate(ST, edges, varargin)

P.smoothing = 2.5e-3;
P = hdsort.util.parseInputs(P, varargin);

if 0
    myhistc = @(ST) histcounts(ST, edges)';
    binnedST_1 = cellfun(myhistc, ST, 'UniformOutput', 0);
    binnedST_2 = cell2mat(binnedST_1);
    
    %binnedST = sum(binnedST_2, 2)';
    binnedST = binnedST_2';
    
    rate = [];
    for ii = 1:size(binnedST,1)
        rate = [rate smooth(binnedST(ii,:), P.smoothing, 'lowess' )];
    end
    
    rate = median(rate')';
    %rate = ~any(~rate')';
    
    %binnedST = binnedST(:)';
    
else
    binnedST = histcounts(cell2mat(ST') , edges);
    rate = smooth(binnedST, P.smoothing, 'lowess' );
end


end