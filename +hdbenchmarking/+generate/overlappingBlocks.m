function [newME, swapPairsEl, blockIdx] = overlappingBlocks(ME)
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

Nmax = 0; newElP = []; newElN = []; swapPairsEl = [];
for b1 = block1
    [newElP_, newElN_, swapPairsEl_] = swapPairs(b1, block2(1));
    if size(swapPairsEl_, 1) > Nmax;
        Nmax = size(swapPairsEl_, 1);
        newElP = newElP_;
        newElN = newElN_;
        swapPairsEl = swapPairsEl_;
    end
end
for b2 = block2
    [newElP_, newElN_, swapPairsEl_] = swapPairs(block1(1), b2);
    if size(swapPairsEl_, 1) > Nmax;
        Nmax = size(swapPairsEl_, 1);
        newElP = newElP_;
        newElN = newElN_;
        swapPairsEl = swapPairsEl_;
    end
end


%%
newME = copy(ME);
newME.electrodePositions = newElP;
newME.electrodeNumbers = newElN;

%% Assign each electrode to a block:
blockIdx = ismember(newME.electrodeNumbers, swapPairsEl(:,1));


    function [newElP, newElN, swapPairsEl] = swapPairs(el1, el2)
        d = elP(el2,:)-elP(el1,:);
        
        newElP = [];
        newElN = [];
        %swapPairsIdx = [];
        swapPairsEl = [];
        for ii = 1:size(elP, 1)
            for jj = 1:size(elP, 1)
                if jj == ii
                    continue;
                end
                if ~any(round(elP(ii,:)-elP(jj,:)+d))
                    newElP = [newElP; elP(ii,:); elP(jj,:)];
                    newElN = [newElN; elN(ii); elN(jj)];
                    swapPairsEl = [swapPairsEl; elN(ii) elN(jj)];
                    %swapPairsIdx = [swapPairsIdx; ii, jj];
                end
            end
        end

    end

end