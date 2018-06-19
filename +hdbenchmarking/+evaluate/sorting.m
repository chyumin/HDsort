function [sortingEvaluation, sortingEvaluationFile] = sorting(...
    artificialUnits, SpikeSortingResult, anaFolder)

%% Prepare output location:
mkdir(anaFolder);
sortingEvaluationFile = fullfile(anaFolder, 'sortingEvaluation.mat');

%%
try
    x
    load(sortingEvaluationFile);
catch
    disp('Analyze the results')
    
    %% Save statistics about each sorting:
    sortingEvaluation.sorted.IDs = SpikeSortingResult.unitIDs;
    sortingEvaluation.gt.IDs = 1:artificialUnits.ST.nCells;
    
    gdf = SpikeSortingResult.getGdf();
    
    %sortingEvaluation.nSpikesInLEG = zeros(numel(LEGs.groups), 1);
    %for gg = 1:numel(LEGs.groups)
    %    sortingEvaluation.nSpikesInLEG(gg) = sum( floor(gdf(:,1)/1000) == gg);
    %end
    
    %% Compare footprints:
    [FPO, fpoFile] = hdbenchmarking.evaluate.footprintOverlap(artificialUnits, SpikeSortingResult, anaFolder);
    
    %% Align spiketrains:
    [alignedSpikeTrains, refractoryPeriod, alignedSpikeTrainFiles] = hdbenchmarking.evaluate.alignSpikeTrains(artificialUnits, SpikeSortingResult, anaFolder);
    
    %% Determint relative values for each unit:
    [sortingEvaluation.matched.TPR, sortingEvaluation.matched.PPV, sortingEvaluation.matched.ERR, sortingEvaluation.matched.FCR, sortingEvaluation.matched.other] ...
        = hdsort.spiketrain.getDetectionInPercents(alignedSpikeTrains);
    sortingEvaluation.matched.meanAMP = [SpikeSortingResult.Units(alignedSpikeTrains.k2f).meanAmplitude];
    
    %% Check whether the footprint overlap results in the same matching
    [m, mi] = min(FPO.shiftedDistances.minDistance);
    sortingEvaluation.matched.sameMatchingWithFPO = mi == alignedSpikeTrains.k2f;
    
    %% Get the amplitude of the inserted unit and set it into relation with the std:
    [mmi, mma, mi_ch, ma_ch] = hdsort.waveforms.tGlobalMinMax(artificialUnits.FP.footprints);
    sortingEvaluation.gt.AMP_abs = abs(mmi)';
    sortingEvaluation.gt.AMP_channels = mi_ch';
    sortingEvaluation.gt.AMP = sortingEvaluation.gt.AMP_abs ./ SpikeSortingResult.noiseStd(sortingEvaluation.gt.AMP_channels)';
    
    %% Get the effective spiking rate of the inserted units:
    totalDuration_s = sum( artificialUnits.ST.lastFrames - artificialUnits.ST.firstFrames + 1) / artificialUnits.ST.samplingRate;
    
    sortingEvaluation.gt.nSpikes = artificialUnits.ST.spikeCount;
    sortingEvaluation.gt.spikingRate_s = sortingEvaluation.gt.nSpikes / totalDuration_s;
    
    %%
    %artificialUnits.maxAMPElectrode = LEGs.electrodeNumbers(mi_ch);
    %artificialUnits.maxAMPLEG = zeros(artificialUnits.nCells, 1);
    %for ku = 1:artificialUnits.nCells
    %    n_ = 0;
    %    for gg = 1:numel(LEGs.groups)
    %        idx = find(LEGs.groups{gg} == artificialUnits.maxAMPElectrode(ku) );
    %        if ~isempty(idx) && sortingEvaluation.nSpikesInLEG(gg)  > n_
    %            n_ = sortingEvaluation.nSpikesInLEG(gg);
    %            artificialUnits.maxAMPLEG(ku) = gg;
    %        end
    %    end
    %end
    %assert(~any(artificialUnits.maxAMPLEG < 1), 'Some LEGs were not correctly detected!')
    
    %sortingEvaluation.percentageSpikesPerLEG = artificialUnits.nSpikes ...
    %    ./ sortingEvaluation.nSpikesInLEG(artificialUnits.maxAMPLEG) * 100;
    
    %if v < 3
    %    assert( ~any(artificialUnits.nSpikes > sortingEvaluation.nSpikesInLEG(artificialUnits.maxAMPLEG) ), ...
    %        'There is a serious problem here!')
    %end
    %%
    unitIDs = SpikeSortingResult.getUnitIDs();
    sortingEvaluation.matched.IDs = unitIDs(alignedSpikeTrains.k2f);
    
    save(sortingEvaluationFile, 'sortingEvaluation')
    
end


end
