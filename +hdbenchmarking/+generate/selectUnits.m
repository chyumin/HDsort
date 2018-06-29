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

%% Select the footprint accorting to the selection criterium
if strcmp(unitselection_parameters.footprintSelectionCriterium, 'random')
    
    %% Choose good footprints based on their 'eligibility' alone:
    originalCellIdxBlock0 = randsample(find(~blockMembership & eligibleUnits), nCellsPerBlock);
    originalCellIdxBlock1 = randsample(find( blockMembership & eligibleUnits), nCellsPerBlock);

elseif strcmp(unitselection_parameters.footprintSelectionCriterium, 'targetamplitude') 
    
    %% Choose good footprints based on their 'eligibility' AND their amplitude:
    idxBlock0 = find(~blockMembership & eligibleUnits);
    idxBlock1 = find( blockMembership & eligibleUnits);
    
    %% Find the amplitudes:
    AMP = [];
    for kk = 1:size(fp, 3);
        m = util.min2D(fp(:,:, kk));
        AMP(kk) = m;
    end
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

assert( sum(ismember(originalCellIdxBlock0, unitselection_parameters.forbiddenUnits)) == 0, 'Gaaaah!')
assert( sum(ismember(originalCellIdxBlock1, unitselection_parameters.forbiddenUnits)) == 0, 'Gaaaah!')

%% Output:
selectedUnits.eligibleUnits = eligibleUnits;
selectedUnits.originalCellIdxBlock0 = originalCellIdxBlock0;
selectedUnits.originalCellIdxBlock1 = originalCellIdxBlock1;
selectedUnits.unitselection_parameters = unitselection_parameters;

end