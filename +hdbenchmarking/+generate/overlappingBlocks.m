function [newME, swapPairsEl, channelIdx, blockIdx] = overlappingBlocks(ME)
% Input:
% DS - Has a MultiElectrode with two high density blocks
%
% Output:
% newDS
% --> newME - MultiElectrode where only those electrodes are retained where the
% two high density blocks overlap (shifted towards each other).
elP = ME.getElectrodePositions;
elN = ME.getElectrodeNumbers;

%% Find the blocks and match the pairs:
idx = kmeans(elP, 2);
block1 = find(idx==1)';
block2 = find(idx==2)';


%% Order blocks:
blockIdx1 = orderBlock(block1);
blockIdx2 = orderBlock(block2);

%% Keep the electrodes that are matching:
combinedBlockIdx = ~isnan( blockIdx1 + blockIdx2);
newBlock1 = blockIdx1(combinedBlockIdx);
newBlock2 = blockIdx2(combinedBlockIdx);

%% Create the new MultiElectrode:
newElP = [elP(newBlock1, :); elP(newBlock2, :)];
newElN = [elN(newBlock1); elN(newBlock2)];

newME = copy(ME);
newME.electrodePositions = newElP;
newME.electrodeNumbers = newElN;

swapPairsEl = [elN(newBlock1)' elN(newBlock2)'];
blockIdx = [zeros(numel(newBlock1), 1); ones(numel(newBlock2), 1)];

for ii = 1:size(newElN,1)
channelIdx(ii) = find( ismember( elN, newElN(ii) ));
end

    function blockIdx = orderBlock(block)
        ux = unique(elP(block, 1));
        uy = unique(elP(block, 2));
        
        %[~, sidx1] = sort( elP(block, 1) );
        %[~, sidx2] = sort( elP(block, 2) );
        
        blockIdx = nan(numel(ux), numel(uy));
        for b = block
            x = find( ismember( ux, elP(b, 1)) );
            y = find( ismember( uy, elP(b, 2)) );
            blockIdx(x,y) = b;
        end
    end

end