function [alignedSpikeTrains, refractoryPeriod, alignedSpikeTrainFiles] = alignSpikeTrains(artificialUnits, SpikeSortingResult, anaFolder)

alignedSpikeTrainFiles = fullfile(anaFolder, ['aligned_spiketrains.mat']);
try
    load(alignedSpikeTrainFiles)
catch
    %%
    disp('Align spiketrains...')
    
    gdf_gt = hdsort.spiketrain.toGdfFrom2DCell(artificialUnits.ST.spikeTrains, SpikeSortingResult.sortingInfo.startTimes(:,1));
    gdf = SpikeSortingResult.getGdf;
    
    maxJitter = 20;
    maxShift = 20;
    maxOverlapDist = 20;
    minimizeError = true;
    
    alignedSpikeTrains = hdsort.spiketrain.alignGDFs(gdf_gt, gdf, ...
        maxJitter, maxShift, maxOverlapDist, ...
        'minimizeError', minimizeError);
    
    %% Check the refractory period violations
    refractoryPeriod.duration_ms = 1.5;
    isih = hdsort.util.isih(gdf, 'binSize_ms', refractoryPeriod.duration_ms);
    refractoryPeriod.violations = isih(: ,1);
    
    %% 
    save(alignedSpikeTrainFiles, 'alignedSpikeTrains', 'refractoryPeriod', '-v7.3')
    
    disp('Done.')
end

end