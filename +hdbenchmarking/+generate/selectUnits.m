function selectedUnits = selectUnits(gdf, footprintT, nCellsPerBlock, blockIdx, unitselection_parameters)

%% Select only units with a spikecount that is larger than the average:
spikeTrains = hdsort.spiketrain.gdf2cell(gdf);
spikeCounts = cellfun( @numel, spikeTrains);
eligibleUnits = spikeCounts > mean(spikeCounts);
eligibleUnits(unitselection_parameters.forbiddenUnits) = false;
eligibleUnits = eligibleUnits(:);

%% Assign an index for the block where the maximum footprint is located:
m1 = max(footprintT, [], 1);
[m11, im11] = max(squeeze(m1), [], 1);
blockMembership = blockIdx(im11);

% blockIdx = [];
% for ii = 1:length(im11)
%     [row, col] = ind2sub(size(swapElectrodePairs),find(swapElectrodePairs==im11(ii)));
%     blockIdx(ii) = col-1;
% end

%% Find the amplitudes:
AMP = [];
for kk = 1:size(footprintT, 3);
    m = hdsort.util.min2D(footprintT(:,:, kk));
    AMP(kk) = m;
end
idxBlock0 = find(~blockMembership & eligibleUnits);
idxBlock1 = find( blockMembership & eligibleUnits);

%% Select the footprint accorting to the selection criterium
if strcmp(unitselection_parameters.footprintSelectionCriterium, 'random')
    
    %% Choose good footprints based on their 'eligibility' alone:
    originalCellIdxBlock0 = randsample(idxBlock0, nCellsPerBlock);
    originalCellIdxBlock1 = randsample(idxBlock1, nCellsPerBlock);

elseif strcmp(unitselection_parameters.footprintSelectionCriterium, 'amplitude_sweep') 
    
    % Block 0:
    amplitudes0 = AMP(idxBlock0);
    [a_, sidx] = sort(amplitudes0);
    block0sidx = sidx(round(linspace(1, length(amplitudes0), nCellsPerBlock)));
    originalCellIdxBlock0 = idxBlock0(block0sidx);
    
    % Block 1:
    amplitudes1 = AMP(idxBlock1);
    [a_, sidx] = sort(amplitudes1);
    block1sidx = sidx(round(linspace(1, length(amplitudes1), nCellsPerBlock)));
    originalCellIdxBlock1 = idxBlock1(block1sidx);
    
elseif strcmp(unitselection_parameters.footprintSelectionCriterium, 'targetamplitude') 
    
    %% Choose good footprints based on their 'eligibility' AND their amplitude:
    targetAMP = mean(AMP)*1.5;
    
    %% Block0:
    [sAMP sIdx] = sort(AMP(idxBlock0));
    [x iMiddle] = min(abs( sAMP - targetAMP) );
    idxRange = (iMiddle-4):(iMiddle+5);
    while any(idxRange == 0)
        idxRange = idxRange + 1;
    end
    originalCellIdxBlock0 = idxBlock0( sIdx(idxRange) ) ;
    
    % Block1:
    [sAMP sIdx] = sort(AMP(idxBlock1));
    [x iMiddle] = min(abs( sAMP - targetAMP) );
    idxRange = (iMiddle-4):(iMiddle+5);
    while any(idxRange == 0)
        idxRange = idxRange + 1;
    end
    originalCellIdxBlock1 = idxBlock1( sIdx(idxRange)) ;
end

%%
assert( numel(originalCellIdxBlock0) == nCellsPerBlock,  'Gaaaah!')
assert( numel(originalCellIdxBlock1) == nCellsPerBlock,  'Gaaaah!')

assert( ~any(ismember(originalCellIdxBlock0, originalCellIdxBlock1)), 'Gaaah!')

assert( sum(ismember(originalCellIdxBlock0, unitselection_parameters.forbiddenUnits)) == 0, 'Gaaaah!')
assert( sum(ismember(originalCellIdxBlock1, unitselection_parameters.forbiddenUnits)) == 0, 'Gaaaah!')

%% Output:
selectedUnits.eligibleUnits = eligibleUnits;
selectedUnits.originalCellIdxBlock0 = originalCellIdxBlock0;
selectedUnits.originalCellIdxBlock1 = originalCellIdxBlock1;
selectedUnits.unitselection_parameters = unitselection_parameters;

end