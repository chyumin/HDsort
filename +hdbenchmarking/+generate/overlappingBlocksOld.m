function [newME, swapPairsEl, blockIdx] = overlappingBlocks(ME, el1, el2, do_plot)
% Input:
% DS - Has a MultiElectrode with two high density blocks
% el1 - corner electrode of first high density block
% el2 - corner electrode of second high density block
%
% Output:
% newDS
% --> newME - MultiElectrode where only those electrodes are retained where the
% two high density blocks overlap (shifted towards each other).
elP = ME.getElectrodePositions;
elN = ME.getElectrodeNumbers;

if nargin < 4
    do_plot = false;
end

d = elP(el2,:)-elP(el1,:);

newElP = [];
newElN = [];
swapPairsIdx = [];
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
            swapPairsIdx = [swapPairsIdx; ii, jj];
        end
    end
end

%newCh = sort(newCh);
if do_plot
    hold on; scatter(newElP(:,1), newElP(:,2))
end
newME = copy(ME);
newME.electrodePositions = newElP;
newME.electrodeNumbers = newElN;

%% Assign each electrode to a block:
blockIdx = ismember(newME.electrodeNumbers, swapPairsEl(:,1));

end