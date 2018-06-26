function [newFootprint] = swapFootprint(footprint, ME, swapPairsEl)
% Swap the waveforms between paired channels
    [nT nCh nU] = size(footprint);
    newFootprint = footprint;
    
    swapPairsCh(:,1) = find(ismember(ME.electrodeNumbers, swapPairsEl(:,1)));
    swapPairsCh(:,2) = find(ismember(ME.electrodeNumbers, swapPairsEl(:,2)));
    
    for u = 1:nU
        newFootprint(:, swapPairsCh(:,1) , u) = footprint(:, swapPairsCh(:,2) , u);
        newFootprint(:, swapPairsCh(:,2) , u) = footprint(:, swapPairsCh(:,1) , u);
    end
end