
%% New file location:
pd = pdefs;
rootFolder = fullfile(pd.mea1kRoot, 'rolandd', 'benchmarking_data_2017');
anaBaseFolder = fullfile(pd.mea1kAnalyzed, 'rolandd', 'benchmarking_data_2017');
originalFolder = fullfile(rootFolder, 'original_recordings');

if v == 1
    datasetName = 'amplitude_sweep01';
elseif v == 2
    datasetName = 'amplitude_sweep02';
elseif v == 3
    datasetName = 'rate_sweep01';
elseif v == 4
    datasetName = 'rate_sweep02';
end

%% Load dataset info file:
for dataSetN = 1:6
    %%
    dataSetNstr = ['dataset' num2str(dataSetN, '%02i')];
    anaFolder = fullfile(anaBaseFolder, datasetName, dataSetNstr);
    recordingsFolder = fullfile(rootFolder, datasetName, dataSetNstr);
    matchedUnitFile = fullfile(anaFolder, 'matchedUnits.mat');
    
    %%
    artificialUnitFile = fullfile(recordingsFolder, 'artificial_units.mat');
    auFileChecksum = Simulink.getFileChecksum(artificialUnitFile);
    
    %%
    try
        uc = load(matchedUnitFile, 'auFileChecksum');
        uc.auFileChecksum
        assert(strcmp(auFileChecksum, uc.auFileChecksum), 'Not the same artificial units!')
    catch
        disp(['Recompute matchedUnitFile for ' datasetName ' ' dataSetNstr])
        
        %%
        mkdir(anaFolder);
        [DH_original, LSA_original, dataSetNstr]   = roland.spikesort_s2.sortOriginal(dataSetN);
        [DH_artificial,  LSA_artificial, dataSetNstr, LEGs, LSA_discarded] = roland.spikesort_s2.sortArtificial(dataSetN, datasetName);
        
        artificialUnits = load(artificialUnitFile);
        
        %% Check AU files:
        for ii = 1:numel(DH_artificial.files.preprocessed)
            auFileChecksum_ = h5read(DH_artificial.files.preprocessed{ii}, '/artificialUnits/checksum');
            auFileChecksum_ = auFileChecksum_{1}(1:end-1);
            assert(strcmp(auFileChecksum_, auFileChecksum), 'Checksum not correct!');
        end
        
        %% Save statistics about each sorting:
        sortingStats.originalIDs = LSA_original.unitIDs;
        sortingStats.artificialIDs = LSA_artificial.unitIDs;
        sortingStats.auIDs = 1:artificialUnits.nCells;
        
        gdf_artificial = LSA_artificial.getGdf;
        sortingStats.nSpikesInLEG = zeros(numel(LEGs.groups), 1);
        for gg = 1:numel(LEGs.groups)
            sortingStats.nSpikesInLEG(gg) = sum( floor(gdf_artificial(:,1)/1000) == gg);
        end
        
        %% Look at the footprint overlaps
        fpoFile = fullfile(anaFolder, ['fpo_new.mat']);
        try
            load(fpoFile)
        catch
            %%
            fp_found = LSA_original.getFootprint();
            vfp_found = mysort.wf.t2v(fp_found);
            FPO.ku = []; FPO.dot_product = [];
            
            fp_found_discarded = LSA_discarded.getFootprint();
            vfp_found_discarded = mysort.wf.t2v(fp_found_discarded);
            FPO_discarded.ku = []; FPO_discarded.dot_product = [];
            for ku = 1:artificialUnits.nCells
                
                fp = artificialUnits.footprints(:,:,ku);
                fp_cut = fp((artificialUnits.cutLeft-20):(artificialUnits.cutLeft+54), :, :);
                vfp = mysort.wf.t2v(fp_cut);
                FPO.ku(:, :, ku) = fp_cut;
                FPO_discarded.ku(:, :, ku) = fp_cut;
                
                %% Dot product:
                FPO.dot_product(ku, :) = vfp * vfp_found';
                FPO_discarded.dot_product(ku, :) = vfp * vfp_found_discarded';
            end
            FPO.found = fp_found;
            FPO_discarded.found = fp_found_discarded;
            
            tic
            [FPO.shiftedDistances, FPO.shiftedDistancesP] = ...
                mysort.wf.tShiftedDistances(FPO.found, FPO.ku, 'maxShift', 10)
            toc
            tic
            [FPO_discarded.shiftedDistances, FPO_discarded.shiftedDistancesP] = ...
                mysort.wf.tShiftedDistances(FPO_discarded.found, FPO_discarded.ku, 'maxShift', 10)
            toc
            save(fpoFile, 'FPO', 'FPO_discarded')
        end
        
        
        if 0
            %% Look at the footprint overlaps
            fpoFile = fullfile(anaFolder, ['fpo.mat']);
            try
                load(fpoFile)
            catch
                %%
                fp_found = LSA_artificial.getFootprint();
                vfp_found = mysort.wf.t2v(fp_found);
                FPO.ku = []; FPO.dot_product = [];
                for ku = 1:artificialUnits.nCells
                    fp = artificialUnits.footprints(:,:,ku);
                    fp_cut = fp((artificialUnits.cutLeft-20):(artificialUnits.cutLeft+54), :, :);
                    
                    vfp = mysort.wf.t2v(fp_cut);
                    
                    FPO.ku(:, :, ku) = fp_cut;
                    FPO.dot_product(ku, :) = vfp * vfp_found';
                end
                FPO.found = fp_found;
                
                tic
                [FPO.shiftedDistances, FPO.shiftedDistancesP] = mysort.wf.tShiftedDistances(FPO.found, FPO.ku, 'maxShift', 10)
                toc
                save(fpoFile, 'FPO')
            end
        end
        
        matchedUnits.FPO = FPO;
        matchedDiscardedUnits.FPO = FPO_discarded;
        clear FPO FPO_discarded
        
        %% Compute and buffer alignment:
        alignedSTFile = fullfile(anaFolder, ['alignedSt.mat']);
        try
            load(alignedSTFile)
            assert(strcmp(auFileChecksum, auFileChecksum_aligned), 'Not the same checksums!')
            disp('Aligned spiketrains loaded')
        catch
            
            %%
            disp('Create alignment...')
            
            gdf_inserted = mysort.spiketrain.toGdfFrom2DCell(...
                artificialUnits.spikeTrains, LSA_artificial.sortingInfo.startTimes);
            gdf_artificial = LSA_artificial.getGdf;
            
            alignedP.maxJitter = 20;
            alignedP.maxShift = 20;
            alignedP.maxOverlapDist = 20;
            alignedP.minimizeError = true;
            
            alignedST = mysort.spiketrain.alignGDFs(gdf_inserted, gdf_artificial, ...
                alignedP.maxJitter, alignedP.maxShift, alignedP.maxOverlapDist, ...
                'minimizeError', alignedP.minimizeError);
            
            auFileChecksum_aligned = auFileChecksum;
            
            %%
            gdf_discarded = DH_artificial.LSA_discarded.getGdf();
            
            alignedST_discarded = mysort.spiketrain.alignGDFs(gdf_inserted, gdf_discarded, ...
                alignedP.maxJitter, alignedP.maxShift, alignedP.maxOverlapDist, ...
                'minimizeError', alignedP.minimizeError);
            
            %%
            ref_period = 1.5;
            isih = util.isih(gdf_artificial, 'binSize_ms', ref_period);
            ref_period_violations = isih(: ,1);
            %%
            save(alignedSTFile, 'alignedST', 'alignedST_discarded', 'alignedP',...
                'auFileChecksum_aligned', 'ref_period', 'ref_period_violations', '-v7.3')
            disp('Alignment finished.')
        end
        
        %% Determint relative values for each unit:
        [matchedUnits.TPR, matchedUnits.PPV, matchedUnits.ER, matchedUnits.FCR, matchedUnits.other] ...
            = mysort.spiketrain.getDetectionInPercents(alignedST);
        matchedUnits.meanAMP = [LSA_artificial.Units(alignedST.k2f).meanAmplitude];
        matchedUnits.totalErrors = alignedST.totErr(:);
        
        [matchedDiscardedUnits.TPR, matchedDiscardedUnits.PPV, matchedDiscardedUnits.ER, matchedDiscardedUnits.FCR, matchedDiscardedUnits.other] ...
            = mysort.spiketrain.getDetectionInPercents(alignedST_discarded);
        matchedDiscardedUnits.meanAMP = [LSA_artificial.Units(alignedST_discarded.k2f).meanAmplitude];
        matchedDiscardedUnits.totalErrors = alignedST_discarded.totErr(:);
        
        %% Footprint similarity & compare the alignment to the discarded
        [m, mi] = min(matchedUnits.FPO.shiftedDistances.minDistance);
        matchedUnits.isMostSimilarFootprintMinDist = mi == alignedST.k2f;
        [m, mi] = max(matchedUnits.FPO.dot_product');
        matchedUnits.isMostSimilarFootprintDotP = mi == alignedST.k2f;
        matchedUnits.isMostSimilarFootprint = matchedUnits.isMostSimilarFootprintMinDist | matchedUnits.isMostSimilarFootprintDotP;
        
        [m, mi_d] = min(matchedDiscardedUnits.FPO.shiftedDistances.minDistance);
        matchedDiscardedUnits.isMostSimilarFootprintMinDist = mi_d == alignedST_discarded.k2f;
        [m, mi_d] = max(matchedDiscardedUnits.FPO.dot_product');
        matchedDiscardedUnits.isMostSimilarFootprintDotP = mi == alignedST_discarded.k2f;
        matchedDiscardedUnits.isMostSimilarFootprint = matchedDiscardedUnits.isMostSimilarFootprintMinDist | matchedDiscardedUnits.isMostSimilarFootprintDotP;
        
        %% Get the amplitude of the inserted unit and set it into relation with the std:
        [mmi, mma, mi_ch, ma_ch] = mysort.wf.tGlobalMinMax(artificialUnits.footprints);
        artificialUnits.AMP_abs = abs(mmi)';
        artificialUnits.AMP_channels = mi_ch';
        artificialUnits.AMP = artificialUnits.AMP_abs ./ LSA_artificial.noiseStd(artificialUnits.AMP_channels)';
        
        % Find LEG in which maximal amplitude was found:
        artificialUnits.nSpikes = sum(artificialUnits.spikeCount, 2);
        
        %%
        artificialUnits.maxAMPElectrode = LEGs.electrodeNumbers(mi_ch);
        artificialUnits.maxAMPLEG = zeros(artificialUnits.nCells, 1);
        for ku = 1:artificialUnits.nCells
            n_ = 0;
            for gg = 1:numel(LEGs.groups)
                idx = find(LEGs.groups{gg} == artificialUnits.maxAMPElectrode(ku) );
                if ~isempty(idx) && sortingStats.nSpikesInLEG(gg)  > n_
                    n_ = sortingStats.nSpikesInLEG(gg);
                    artificialUnits.maxAMPLEG(ku) = gg;
                end
            end
        end
        assert(~any(artificialUnits.maxAMPLEG < 1), 'Some LEGs were not correctly detected!')
        
        matchedUnits.percentageSpikesPerLEG = artificialUnits.nSpikes ...
            ./ sortingStats.nSpikesInLEG(artificialUnits.maxAMPLEG) * 100;
        
        if v < 3 
            assert( ~any(artificialUnits.nSpikes > sortingStats.nSpikesInLEG(artificialUnits.maxAMPLEG) ), ...
                'There is a serious problem here!')
        end
        %%
        IDs = LSA_artificial.unitIDs(alignedST.k2f)
        artificialUnits.detectedLEG = floor(IDs/1000);
        
        assert(~any(artificialUnits.nSpikes > sortingStats.nSpikesInLEG(artificialUnits.detectedLEG)), 'There are more inserted spikes than found spikes!')
        
        matchedUnits.percentageSpikesPerDetectedLEG = artificialUnits.nSpikes ...
            ./ sortingStats.nSpikesInLEG(artificialUnits.detectedLEG) * 100;
        %%
        
        
        save(matchedUnitFile, 'sortingStats', 'matchedUnits', 'matchedDiscardedUnits', ...
            'artificialUnits', 'auFileChecksum')
        
        %%
        clear sortingStats matchedUnits matchedDiscardedUnits artificialUnits auFileChecksum
        
    end
end

%% ------------------------------------------------------------------------
if 0
    %%
    v = 3; dataSetN = 1;
    roland.spikesort_s2.load_results
    
    fullSTcomparisonFile = fullfile(anaFolder, ['fullSTcomparison.mat']);
    
    try
        load(fullSTcomparisonFile)
    catch
        gdf_inserted = mysort.spiketrain.toGdfFrom2DCell(artificialUnits.spikeTrains, LSA_artificial.sortingInfo.startTimes);
        gdf_artificial = LSA_artificial.getGdf;
        gdf_original = LSA_original.getGdf;
        
        %%
        alignedP.maxJitter = 20;
        alignedP.maxShift = 20;
        alignedP.maxOverlapDist = 20;
        alignedP.minimizeError = true;
        
        R = mysort.spiketrain.alignGDFs(gdf_original, gdf_artificial, ...
            alignedP.maxJitter, alignedP.maxShift, alignedP.maxOverlapDist, ...
            'minimizeError', alignedP.minimizeError);
        
        save(fullSTcomparisonFile, 'R', 'alignedP')
    end
    
    %%
    %     maxJitter = 0;
    %     maxShift = 0;
    %
    %     gdf_original_x   = gdf_original;   gdf_original_x(:,1) = 1;
    %     gdf_artificial_x = gdf_artificial; gdf_artificial_x(:,1) = 1;
    %
    %     st_original = {gdf_original(:,2)};
    %     st_artificial = {gdf_artificial(:,2)};
    %
    %     R = mysort.spiketrain.compare(st_original, st_artificial, maxJitter, maxShift, false)
    
end
