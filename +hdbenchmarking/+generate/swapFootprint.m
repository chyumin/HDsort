function [newFootprint, blockIdx] = swapFootprint(footprint, ME, swapPairsEl)
% Swap the waveforms between paired channels
    [nT nCh nU] = size(footprint);
    newFootprint = footprint;
    for ii = 1:size(swapPairsEl,1)
        swapPairsEl_(ii,1) = find(ME.electrodeNumbers(:) == swapPairsEl(ii,1));
        swapPairsEl_(ii,2) = find(ME.electrodeNumbers(:) == swapPairsEl(ii,2));
    end
    swapPairsEl = swapPairsEl_;
    for u = 1:nU
        newFootprint(:, swapPairsEl(:,1) , u) = footprint(:, swapPairsEl(:,2) , u);
        newFootprint(:, swapPairsEl(:,2) , u) = footprint(:, swapPairsEl(:,1) , u);
    end
    
    %% Assign an index for the block where the maximum footprint is located:
    m1 = max(newFootprint, [], 1);
    [m11 im11] = max(squeeze(m1), [], 1);
    blockIdx = [];
    for ii = 1:length(im11)
        [row, col] = ind2sub(size(swapPairsEl),find(swapPairsEl==im11(ii)));
        blockIdx(ii) = col-1;
    end
end