function [S P] = sortScript(DS, dpath, name, varargin)
    S.STOP_ME_BECAUSE_I_AM_SLOW = false;
    % Artefact Detection
    P.artefactDetection.use = 0;
    P.artefactDetection.width = 350;
    P.artefactDetection.threshold = size(DS,2)*41; % should be sampling frequency dependent!
    
    % Spike Detection
    P.spikeDetection.method = '-';
    P.spikeDetection.thr = 3.5;
    P.spikeDetection.minDist = 25;
    P.spikeDetection.maxDataLength = [];
    P.spikeDetection.mergeSpikesMaxDist = 18; % do not go below 13 for HDMEA!
    P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan = 1500000;
    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude = [0 0 0]; % should remove atrifacts that are oscillations (successive detections of peaks with alternating signs). First parameter is N, second the minHeight, third the maxDist within a group
    
    % Spike cutting
    %P.spikeCutting.maxSpikes = 100000;
    % Do not use (i.e. make it smaller) this parameter unless you know what
    % you are doing. This will have the consequence that not all spikes
    % will be cut and thus also not BOTM-matched. This will cause troubles
    % for the HD postprocessing "ana.processLocalSortings". In this case
    % either the full BOTM needs to be run or a batched version of the cut-
    % and-match function implemented.
    P.spikeCutting.maxSpikes = 200000000000;
    
    P.spikeCutting.Tf = 71;
    P.spikeCutting.cutLeft = 20;
    P.spikeCutting.chunkSize = 10000;
    P.spikeCutting.blockwise = true;
    
    % Spike Alignment
    P.spikeAlignment.method = 'onUpsampledMean';
    P.spikeAlignment.maxSpikes = 50000;
    P.spikeAlignment.Tf = 35;
    P.spikeAlignment.cutLeft = 12;
    P.spikeAlignment.initAlignment = []; %'-' % Do not use init alignment
    P.spikeAlignment.maxIdx = P.spikeAlignment.cutLeft + 1;
    
    P.spikeAlignment.maxIterations = 30;
    
    % Noise estimation    
    P.noiseEstimation.maxSamples = 1000000;
    P.noiseEstimation.maxSamplesPerEpoch = 500000;
    P.noiseEstimation.minDistFromSpikes = 60;
    
    % Feature extraction
    P.featureExtraction.Tf = 20;
    P.featureExtraction.cutLeft = 8;
    P.featureExtraction.nDims = 6;
    
    % Clustering
    P.clustering.maxSpikes = 30000;
    P.clustering.meanShiftBandWidth = []; % default is sqrt( P.clustering.meanShiftBandWidth * P.featureExtraction.nDims )
    P.clustering.meanShiftBandWidthFactor = 1.1;
%     P.clustering.bandwidthMultMaxDistSpeedMS = 5;
    P.clustering.minSpikesPerCluster = 10;
    
    % Template Matching on Cut Spikes
    P.templateMatchingCut.prior = .0001;
    P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];
    
    % Merge Templates after template matching
    P.mergeTemplates.merge = 1;
	P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .96;
%     P.mergeTemplates.ifMeanDistSmaller = .1;    % in units std of noise
    P.mergeTemplates.ifMaxDistSmaller = 2.5;      % in units std of noise
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 25;

    % Run BOTM Template Matching on whole data
    P.botm.run = 0;
    P.botm.Tf = 55;
    P.botm.cutLeft = 10;
    P.botm.prior = .0001;

    P = hdsort.util.parseInputs(P, varargin, 'error');

    % Checks
    assert(P.botm.run==0 ||P.botm.cutLeft <= P.spikeCutting.cutLeft, 'Spike Cut CutLeft must be greater or equal to Botm cut!');
    assert(P.botm.run==0 || P.botm.Tf <= P.spikeCutting.Tf, 'Spike Cut Tf must be greater or equal to Botm Tf!');
    assert(P.botm.run==0 ||P.spikeCutting.cutLeft - P.botm.cutLeft + P.botm.Tf <= P.spikeCutting.Tf, 'Spike Cut for Botm out of bounds!');
    assert(P.spikeCutting.cutLeft >= P.spikeAlignment.cutLeft, 'Cannot align further left than cut !');
    assert(P.featureExtraction.Tf < P.spikeAlignment.Tf, 'Feature extraction spikeforms must be shorter than the alignment spikeforms!')
    assert(P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft >= 0, 'Connot prewhiten spikes with this cutleft!');
    % Save header
    readme = 'The files in this folder were created by the sort.m script.';
    
    % Define global names
    S.P = P;
    S.name = name;
    S.srate = DS.getSamplesPerSecond();
    S.nC = size(DS,2);
    S.dpath = dpath;
    S.dprefix = fullfile(dpath, name);
    
    if ~exist(S.dpath, 'file')
        fprintf('Output directory does not exists. Trying to create...')
        mkdir(S.dpath);
        fprintf(' success.\n')
    end    
    
    % Define file names
    S.files.artefact_file          = [S.dprefix '.010artefacts.mat'];
    S.files.spike_det_file         = [S.dprefix '.020spikes_det.mat'];
    S.files.spike_det_merged_file  = [S.dprefix '.030spikes_det_merged.mat'];
    S.files.spike_cut_file         = [S.dprefix '.040spikes_cut.mat'];
    S.files.spike_aligned_file     = [S.dprefix '.050spikes_aligned.mat'];
    S.files.cov_file               = [S.dprefix '.060cov.mat'];
    S.files.prewh_spike_file       = [S.dprefix '.070spikes_prewhitened.mat'];
    S.files.fet_spike_file         = [S.dprefix '.080spikes_features.mat'];
    S.files.meanshift_spike_file   = [S.dprefix '.090clusters_meanshift.mat'];
    S.files.botm_matching_file     = [S.dprefix '.100botm_matching.mat'];    
    S.files.merge_ms_clusters      = [S.dprefix '.110clusters_meanshift_merged.mat'];     
    S.files.botm_file              = [S.dprefix '.120botm.mat'];
    S.files.botm_templates_file    = [S.dprefix '.130botm_templates.mat'];
    S.files.resArtefact_file       = [S.dprefix '.140residual_artefacts.mat'];
    S.files.botm_templates_cleaned_file = [S.dprefix '.150botm_cleaned_templates.mat'];
    save([fullfile(dpath, name) '.P.mat'], 'S', 'readme');
    
    % Define constants
    nC = size(DS,2);    

    % RUN
    runArtefactDetection(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runSpikeDetection(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    mergeDetectedSpikes(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    noiseEstimation(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    cutSpikes(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    alignSpikes();     if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    prewhitenSpikes(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    fetExtraction(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runMeanShift(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runBOTMMatching(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    mergeClustersMS(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runBOTM(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
%     computeBOTMTemplates(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
%     residualArtefactDetection();  if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end   
%     computeBOTMCleanedTemplates(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
%     mergeClustersBotm(); if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
   
    disp('############################');
    disp('All done.!');
    
    
    %----------------------------------------------------------------------
    % Detected sp:  +++++ ++ +++++++++++++ +++++++ ++++++++++ ++++++ +++++
    % Cut spikes :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    % Aligned sp :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Prewhitend :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Features   :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Clustered  :     +     +           +     + +                 + +
    % Matched    :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    % BOTM       :  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    %% --------------------------------------------------------------------
    function runArtefactDetection()
        artefactDetection.epochs = [];        
        if P.artefactDetection.use
            try
                load(S.files.artefact_file, 'artefactDetection');
                disp('artefacts already detected');
            catch
                disp('detecting artefacts...');
                A = ana.moritzheimdahl.ArtefactDetector(DS, P.artefactDetection.width, P.artefactDetection.threshold);
                artefactDetection.epochs = A.getEpochs();
                saveIfDoesNotExist(S.files.artefact_file);                 
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
                save(S.files.artefact_file, 'artefactDetection');
                disp('Done.')    
            end
        end
        S.artefactDetection = artefactDetection; clear artefactDetection;
    end

    %% --------------------------------------------------------------------
    function runSpikeDetection()
        try
            spikeDetection = [];            
            load(S.files.spike_det_file);
            disp('spikes already detected');            
        catch
            disp('Detecting spikes...');
            L = P.spikeDetection.maxDataLength;
            spikeDetection.spikesDetectedUp = {};
            spikeDetection.pks_up = {};
            if ~isempty(strfind(P.spikeDetection.method, '+'))
                [spikeDetection.spikesDetectedUp spikeDetection.pks_up] = ...
                    DS.detectSpikes('Len', L, 'energyfun', @(x) x, 'minPeakDistance', P.spikeDetection.minDist, 'thr', P.spikeDetection.thr);
            end        
            spikeDetection.spikesDetectedDown = {};
            spikeDetection.pks_down = {};
            if ~isempty(strfind(P.spikeDetection.method, '-'))
                [spikeDetection.spikesDetectedDown spikeDetection.pks_down] = ...
                    DS.detectSpikes('Len', L, 'energyfun', @(x) -x, 'minPeakDistance', P.spikeDetection.minDist, 'thr', P.spikeDetection.thr);
            end
            saveIfDoesNotExist(S.files.spike_det_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.spike_det_file, 'spikeDetection');
            disp('Done.')    
        end
        S.spikeDetection = spikeDetection; clear spikeDetection;
    end

    %% --------------------------------------------------------------------
    function mergeDetectedSpikes()
        try
            spikeDetectionMerged = [];
            load(S.files.spike_det_merged_file);
            disp('spikes already merged');
        catch
            disp('Merging spikes trains...');
            allspikes = [cell2mat(S.spikeDetection.spikesDetectedUp)   cell2mat(S.spikeDetection.pks_up);
                         cell2mat(S.spikeDetection.spikesDetectedDown) cell2mat(S.spikeDetection.pks_down)];
            nSp = size(allspikes,1);
            if nSp == 0
                disp('No Spikes found in this group!!!!')
                spikeDetectionMerged.ts = [];
                spikeDetectionMerged.allspikes = [];
                S.spikeDetectionMerged = spikeDetectionMerged; 
                fprintf('Found 0 spikes after merging.\n');
                save(S.files.spike_det_merged_file, 'spikeDetectionMerged');
                return
            end
            % add the channel information to allspikes
            allspikes(nSp,3) = 0;
            nextIdx = 1;
            LL = max(length(S.spikeDetection.spikesDetectedUp),...
                     length(S.spikeDetection.spikesDetectedDown));
            for i=1:LL
                nSpUp = 0;
                nSpDo = 0;
                if ~isempty(S.spikeDetection.spikesDetectedUp)
                    nSpUp = length(S.spikeDetection.spikesDetectedUp{i});
                end
                if ~isempty(S.spikeDetection.spikesDetectedDown)
                    nSpDo = length(S.spikeDetection.spikesDetectedDown{i});
                end
                allspikes(nextIdx:nextIdx+nSpUp+nSpDo-1,3) = i;
                nextIdx = nextIdx+nSpUp+nSpDo;
            end
            allspikes = sortrows(allspikes,1);
            fprintf('%d spikes found total before merging\n', nSp);
            % Check for oscillation groups and remove those before
            % artifact removal
            if P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1) > 0
                allspikes = spiketrain.findOscillationGroups(allspikes, ...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(3), ...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1),...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(2));
            end
            % remove spikes in artefact epochs
            removeIdx = hdsort.epoch.findPointsInEpochs(allspikes(:,1), S.artefactDetection.epochs)>0;               
            allspikes(removeIdx,:) = [];
            allspikes(abs(allspikes(:,2))>P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan,:) = [];
            fprintf('%d spikes after artefact removal (%d removed)\n', size(allspikes,1), nSp-size(allspikes,1));

            allspikes  = hdsort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, P.spikeDetection.mergeSpikesMaxDist);
            
            % remove spikes that are too close to the beginning of the file
            removeTooEarly = find(allspikes(:,1) < P.spikeCutting.Tf);
            if ~isempty(removeTooEarly)
                allspikes(removeTooEarly,:) = [];
            end
            % remove spikes that are too close to the end of the file
            removeTooLate = find(allspikes(:,1) > size(DS,1) - P.spikeCutting.Tf);
            if ~isempty(removeTooLate)
                allspikes(removeTooLate,:) = [];
            end
            
            spikeDetectionMerged.allspikes = allspikes;
            spikeDetectionMerged.ts = allspikes(:, 1);
            saveIfDoesNotExist(S.files.spike_det_merged_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.spike_det_merged_file, 'spikeDetectionMerged');
            disp('Done.')    
        end
        S.spikeDetectionMerged = spikeDetectionMerged; 
        fprintf('Found %d spikes after merging.\n', length(spikeDetectionMerged.ts));
    end

    %% --------------------------------------------------------------------
    function cutSpikes()
        fname = [S.files.spike_cut_file '_cutSpikes.bin'];
        try
            spikeCut = [];
            if exist(S.files.spike_cut_file, 'file')
                
                % Check if file is complete (if the variable does not
                % exist, this might be an old file and we simply load it)
                if hdsort.file.hdf5.exist(S.files.spike_cut_file, '/bFileIsComplete')
                    bFileIsComplete = hdf5read(S.files.spike_cut_file, '/bFileIsComplete');
                    if ~bFileIsComplete
                        warning('Spike cutting was not complete, overwriting old files');
                        delete(S.files.spike_cut_file);
                        try
                            delete(fname);
                        catch
                        end
                    end
                else
                    warning('Old spike cutting file, deleting...');
                    delete(S.files.spike_cut_file);
                    try
                        delete(fname);
                    catch
                    end
                end
            end
            load(S.files.spike_cut_file);
            spikeCut.wfs = hdsort.file.hdf5.matrix(fname, '/wfs', 1);
            
            if ~any(any(spikeCut.wfs(:,:) ~= 0.0))
                try
                    delete(S.files.spike_cut_file)
                catch
                end
                try
                    delete(fname)
                catch
                end
                try
                    delete(S.files.spike_aligned_file)
                    delete(S.files.prewh_spike_file)
                    delete(S.files.fet_spike_file)
                    delete(S.files.meanshift_spike_file)
                    delete(S.files.botm_matching_file)
                    delete(S.files.merge_ms_clusters)
                catch
                    disp('Not all files could be deleted!')
                end
                error('All traces are zero, recut them!')
            end
            
            disp('spikes already cut');
        catch
            disp('Cutting spikes...');
            nSpikesDetected = length(S.spikeDetectionMerged.ts);
            if nSpikesDetected > P.spikeCutting.maxSpikes
                cutIdx = randperm(nSpikesDetected);
                spikeCut.cutIdx = sort(cutIdx(1:P.spikeCutting.maxSpikes));
            else
                spikeCut.cutIdx = 1:nSpikesDetected;
            end
            if ~isempty(spikeCut.cutIdx)
                cut_ts = S.spikeDetectionMerged.ts(spikeCut.cutIdx);
                nSP = length(cut_ts);
                nC = size(DS,2);
                L = nC*P.spikeCutting.Tf;
                dims = [nSP L];
                maxDims = dims;
                h5type = 'H5T_NATIVE_FLOAT';
                chunkDims = [];
                deflation = 0;
                
                %% CHECK ALREADY IF SPIKES NEED TO BE ALIGNED, IF SO, KEEP IN MEMORY
                
                if P.spikeAlignment.maxSpikes < nSP
                    alignIdx = randperm(nSP);
                    spikeCut.alignIdx = sort(alignIdx(1:P.spikeAlignment.maxSpikes));
                else
                    spikeCut.alignIdx = 1:nSP;
                end
                spikeCut.nSpikesToAlign = length(spikeCut.alignIdx);
                spikeCut.unalignedwfs = zeros(spikeCut.nSpikesToAlign, L);
                
                spikeCut.wfs = hdsort.file.hdf5.matrix(fname, '/wfs', 0, dims, maxDims, h5type, chunkDims, deflation);

                chk = hdsort.util.Chunker(nSP, 'chunkSize', P.spikeCutting.chunkSize, ...
                                 'chunkOverlap', 0);
                while chk.hasNextChunk()
                    fprintf('Cutting Spikes Chunk %d of %d\n', chk.currentChunk, chk.nChunks);
                    [chunkOverlap, chunk, chunkLen] = chk.getNextChunk();
                    chidx = chunk(1):chunk(2);
                    wfs = single(DS.getWaveform(cut_ts(chidx),...
                            P.spikeCutting.cutLeft, P.spikeCutting.Tf)); 
                    spikeCut.wfs(chidx,1:L) = wfs;
                        
                    % Store some of the waveforms in memory for later
                    % alignment
                    spikes2AlignInThisChunk = spikeCut.alignIdx >= chunk(1) & spikeCut.alignIdx <= chunk(2);
                    idxIntoWfs = spikeCut.alignIdx(spikes2AlignInThisChunk);
                    spikeCut.unalignedwfs(spikes2AlignInThisChunk,:) = double(wfs(idxIntoWfs-chunk(1)+1,:));
                end
                
            else
                spikeCut.wfs = [];
            end
            % Create File and complete variable
            bFileIsComplete = int32(0);
            saveIfDoesNotExist(S.files.spike_cut_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.spike_cut_file, 'bFileIsComplete', '-v7.3');
            save(S.files.spike_cut_file, 'spikeCut', '-append');
            % Flag File as complete
            bFileIsComplete = int32(1);
            save(S.files.spike_cut_file, 'bFileIsComplete', '-append');
            disp('Done.')    
        end
        S.spikeCut = spikeCut; clear spikeCut;
    end

    %% --------------------------------------------------------------------
    function noiseEstimation()
        try
            noise = [];
            load(S.files.cov_file);
            disp('Cov already estimated...');
        catch
            disp('Estimating Cov...');
            s1 = S.spikeDetectionMerged.ts-P.noiseEstimation.minDistFromSpikes;
            s2 = S.spikeDetectionMerged.ts+P.noiseEstimation.minDistFromSpikes;
            if isempty(P.spikeDetection.maxDataLength)
                L = size(DS,1);
            else
                L  = P.spikeDetection.maxDataLength;
            end
            noise.epochs = hdsort.epoch.flip(hdsort.epoch.merge([s1 s2]), L);

            maxTf = max([P.botm.Tf P.featureExtraction.Tf P.spikeAlignment.Tf]); %P.spikeCutting.Tf
            Cest = hdsort.noise.Covest2(DS, 'maxLag', maxTf,...
                'maxSamplesPerEpoch', P.noiseEstimation.maxSamplesPerEpoch, ...
                'maxSamples', P.noiseEstimation.maxSamples, ...
                'noiseEpochs', noise.epochs, 'forceMethod', 'xcorr');
            noise.C_time_cut = hdsort.noise.ccol2Cte(Cest.CCol, maxTf);
            noise.C_time_aligned = hdsort.noise.ccol2Cte(Cest.CCol, P.spikeAlignment.Tf);
            noise.meanNoiseStd = sqrt(mean(diag(noise.C_time_cut)));
            noise.CestS = Cest.toStruct();    
            saveIfDoesNotExist(S.files.cov_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.cov_file, 'noise');
            disp('Done.')    
        end
        S.noise = noise; clear noise;
    end

    %% --------------------------------------------------------------------
    function alignSpikes()
        % Make sure to align the spike maximall in the window of the dead
        % time of the spike detector. Otherwise, we detect to spikes which
        % are close to each other and then we shift them to effectively
        % detect the same spike twice.
        if ~isempty(P.spikeAlignment.method); 
            try
                spikeAligned = [];
                load(S.files.spike_aligned_file);
                disp('Spikes already aligned');
            catch
                disp('Aligning spikes...');
                idxstart = 1 + P.spikeCutting.cutLeft - P.spikeAlignment.cutLeft;
                idxstop  = idxstart + P.spikeAlignment.Tf - 1;
                nSpCut = size(S.spikeCut.wfs,1);
                
                %% %%%%%%%%%%%%%%%%%
                if nSpCut == 0
                    spikeAligned.tau = [];
                    spikeAligned.maxIdx = [];
                    spikeAligned.restrictToIdx =[];
                    spikeAligned.wfs = [];
                    spikeAligned.alignIdx = [];
                else
                    spikeAligned.alignIdx = S.spikeCut.alignIdx;
                    spikeAligned.unalignedwfs = double(hdsort.waveforms.vSubsel(S.spikeCut.unalignedwfs,...
                        nC, idxstart:idxstop));    

                    spikeAligned.tau = zeros(size(spikeAligned.unalignedwfs,1),1);
                    spikeAligned.maxIdx = P.spikeAlignment.cutLeft;
                    spikeAligned.idxstart = idxstart;
                    spikeAligned.idxstop = idxstop;

                    spikeAligned.restrictToIdx = spikeAligned.maxIdx-2:spikeAligned.maxIdx+5;
                    tic
                    if strcmp(P.spikeAlignment.method, 'onMax')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnMax(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMean')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnUpsampleMean(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', P.spikeAlignment.maxIdx,...
                            'maxShiftPerIter', 3,...
                            'maxIter', P.spikeAlignment.maxIterations, ...
                            'initAlignment', P.spikeAlignment.initAlignment);   
%                             'restrictToIdx',  spikeAligned.restrictToIdx,...                        
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.alignWaveformsUpsampleMax(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.alignWaveformsUpsampleMax(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2); 
                        spikeAligned.wfs = -spikeAligned.unalignedwfs;
                    elseif strcmp(P.spikeAlignment.method, 'onMin')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnMax(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(P.spikeAlignment.method, 'onAverageMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.vAlignOnAverageMaxSample(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(P.spikeAlignment.method, 'onAverageMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.vAlignOnAverageMaxSample(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(P.spikeAlignment.method, 'none')
                        spikeAligned.wfs = spikeAligned.unalignedwfs;
                    else
                        error(['unkown alignement method' P.spikeAlignment.method]);
                    end
                    toc
                end
%                  figure, plot(spikeAligned.wfs(1:1000,:)')
                % cut away alignement artefacts at the ends
%                 spikeAligned.unalignedwfs([1:max(1,max(spikeAligned.tau)) end+min(0,min(spikeAligned.tau))+1:end],:) = [];
                saveIfDoesNotExist(S.files.spike_aligned_file);                 
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
                save(S.files.spike_aligned_file, 'spikeAligned', '-v7.3');
                disp('Done.')
            end
        end
        S.spikeAligned = spikeAligned; clear spikeAligned;        
    end

    %% --------------------------------------------------------------------
    function prewhitenSpikes()
        try
            spikePrewhitened = [];
            load(S.files.prewh_spike_file);
            disp('Spikes already prewhitened...');
        catch
            disp('Prewhitening spikes...');
            % Cut the aligned waveforms to get rid of alignement artefacts
            % and reduce to the final waveform that will be used for
            % feature extraction
                idxstart = 1 + P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft;
                idxstop  = idxstart + P.featureExtraction.Tf - 1;   
                ccol_loaded = S.noise.CestS.CCol;
                ccol_loaded(1:nC, 1:nC) = ccol_loaded(1:nC, 1:nC) + eye(nC) * .1 * mean(diag(ccol_loaded(1:nC, 1:nC)));                
                spikePrewhitened.ccol_loaded = ccol_loaded/2;
                spikePrewhitened.C = hdsort.noise.ccol2Cte(spikePrewhitened.ccol_loaded, P.featureExtraction.Tf);
                spikePrewhitened.U = chol(spikePrewhitened.C);            
            if isempty(S.spikeAligned.wfs)
                % Prewhiten
                spikePrewhitened.wfs = [];
            else
                spikePrewhitened.wfs = hdsort.waveforms.vSubsel(S.spikeAligned.wfs, nC, idxstart:idxstop);      
                % Build the noise covariance matrix and load it
                spikePrewhitened.wfs = spikePrewhitened.wfs/spikePrewhitened.U; 
            end
            saveIfDoesNotExist(S.files.prewh_spike_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.prewh_spike_file, 'spikePrewhitened', '-v7.3');
            disp('Done.')    
        end
        S.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
    end

    %% --------------------------------------------------------------------
    function fetExtraction()
        try
            spikeFeatures = [];
            load(S.files.fet_spike_file );
            disp('Features already calculated...');
        catch
            disp('Calculating features...');
            if ~isempty(S.spikePrewhitened.wfs)
                [spikeFeatures.X, spikeFeatures.PCs] = hdsort.util.dimReductionPCA(S.spikePrewhitened.wfs,...
                    P.featureExtraction.nDims, [], 3*1000000);
            else
                spikeFeatures.X = [];
            end
            saveIfDoesNotExist(S.files.fet_spike_file);
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.fet_spike_file, 'spikeFeatures', '-v7.3');
            disp('Done.')    
        end
        S.spikeFeatures = spikeFeatures; clear spikeFeatures;
    end
    
    %% --------------------------------------------------------------------
    function runMeanShift()
        try
            clustering = [];
            load(S.files.meanshift_spike_file);
            disp('Already clustered...');
        catch
            disp('Mean Shift Clustering...');
            if ~isempty(P.clustering.meanShiftBandWidth)
                clustering.bandwidth = P.clustering.meanShiftBandWidth;
            else
                clustering.bandwidth = sqrt(P.clustering.meanShiftBandWidthFactor*P.featureExtraction.nDims);
            end
            nSpFet = size(S.spikeFeatures.X,1);
            S.maxNBandwidthIncreases = 1;
            S.bandwidthIncreaseFactor = 1.3;
            if ~isempty(S.spikeFeatures)
                % Select a random subset of spikes:
                if P.clustering.maxSpikes < nSpFet
                    clusterIdx = randperm(nSpFet);
                    clustering.clusterIdx = sort(clusterIdx(1:P.clustering.maxSpikes));
                else
                    clustering.clusterIdx = 1:nSpFet;
                end
                X = S.spikeFeatures.X(clustering.clusterIdx,:)';
                
                t_ = tic;
    %             [clustCent,point2cluster,clustMembsCell] = meanshift.MeanShiftCluster(X, clustering.bandwidth);
                [clustCent,point2cluster,clustMembsCell] = meanshift.MeanShiftClusterIncreaseBW(X, clustering.bandwidth, 0 , P.clustering.minSpikesPerCluster, S.maxNBandwidthIncreases, S.bandwidthIncreaseFactor);
                t_ = toc(t_);
                fprintf('Mean Shift Clustering took %.1f sec\n', t_);
                [clustering.ids clustering.clusterCenter] = meanshift.MeanShiftClusterBundleResult(X', clustMembsCell, P.clustering.minSpikesPerCluster);
                %*** clustering.ids==0 --> noise cluster ***
                clustering.classes = unique(clustering.ids);
                clusteredAlignedSpikes = S.spikeAligned.wfs(clustering.clusterIdx,:);
                clusteredCutSpikes = S.spikeCut.wfs(S.spikeAligned.alignIdx(clustering.clusterIdx),:);
                
                clustering.templatesCut = hdsort.util.calculateClassMeans(clusteredCutSpikes, clustering.ids, 'usemedian'); % templates in order of unique(clustering.ids)
                clustering.templatesAligned = hdsort.util.calculateClassMeans(clusteredAlignedSpikes, clustering.ids, 'usemedian');
            else
                clustering.ids = [];
                clustering.clusterCenter = [];
                clustering.classes = [];
                clustering.templatesCut  = [];
                clustering.templatesAligned  = [];               
            end
            saveIfDoesNotExist(S.files.meanshift_spike_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.meanshift_spike_file, 'clustering', '-v7.3');
        end
        S.clustering = clustering; clear clustering;
    end
    %% --------------------------------------------------------------------
    function runBOTMMatching()
        recalc = false;
        try
            clusteringMatched = [];
            load(S.files.botm_matching_file);
            disp('BOTM matched with lda...');
        catch
            recalc = true;
        end
        
        if recalc
            disp('BOTM matching...');
            nSpCut = size(S.spikeCut.wfs,1);
            if nSpCut == 0
                clusteringMatched.ids = [];
                clusteringMatched.template_ids = [];
                clusteringMatched.spikeCutAligned = [];
                clusteringMatched.templates = [];
                clusteringMatched.maxTausPerSpike = [];
                clusteringMatched.ts = [];
                clusteringMatched.maxTausPerSpikeAndFilter  = [];
                clusteringMatched.tauRange = [];
                clusteringMatched.subTauRange  = [];
                clusteringMatched.tauIdx = [];
                clusteringMatched.subTauIdx  = [];
                clusteringMatched.TupShiftDown = [];
                clusteringMatched.FupShiftDown  = [];
                clusteringMatched.maxTausPerSpike = [];
                clusteringMatched.maxTausIdxPerSpike = [];
                clusteringMatched.maxSubTausIdxPerSpike = [];
                clusteringMatched.ids = [];
                clusteringMatched.template_ids = [];
                clusteringMatched.templatesIDs  = [];
                clusteringMatched.ampsChIdx = [];
                clusteringMatched.amps = [];
            else
                clusteringMatched.ids = zeros(nSpCut,1);
                clusteringMatched.template_ids = 0;
                clusteringMatched.spikeCutAligned = [];
                clusteringMatched.templates = [];
                clusteringMatched.maxTausPerSpike = zeros(nSpCut,1);
                clusteringMatched.ts = S.spikeDetectionMerged.ts(S.spikeCut.cutIdx);

                if size(S.clustering.templatesCut,1) < 2
                    disp('Cannot run BOTM Matching, not enough templates!');
                else
                    T = S.clustering.templatesAligned(2:end,:); % remove the noise template
                    nT = size(T,1);
                    C = (S.noise.C_time_aligned +diag(diag(S.noise.C_time_aligned)))/2;

                    % Do matching only on temporal window defined by the aligned
                    % templates
                    idxstart = 1 + P.spikeCutting.cutLeft - P.spikeAlignment.cutLeft;
    %                 idxstop  = idxstart + P.spikeAlignment.Tf - 1;

                    [D maxTaus TupShiftDown FupShiftDown tauRange tauIdx subTauRange subTauIdx] =hdsort.waveforms.vTemplateMatching(...
                        S.spikeCut.wfs, T, nC, idxstart, 'maxShift', 5, 'upsample', 5, 'noiseCovariance', C, 'chunkSize', P.spikeCutting.chunkSize);
                    assert(any(D(:)~=0.0), 'SpikeCut file corrupt!')
                    F = FupShiftDown(:,:,1);
                    EF = diag(T*F');                     % compute energies
                    Prior  = P.templateMatchingCut.prior;                
                    DISCR = D - .5 * repmat(EF', nSpCut, 1)  + log(Prior);  % compute botm dicriminant
                    clusteringMatched.maxTausPerSpikeAndFilter = maxTaus;
                    clusteringMatched.tauRange = tauRange;
                    clusteringMatched.subTauRange = subTauRange;
                    clusteringMatched.tauIdx = tauIdx;
                    clusteringMatched.subTauIdx = subTauIdx;
                    clusteringMatched.TupShiftDown = TupShiftDown;
                    clusteringMatched.FupShiftDown = FupShiftDown;

                    [MaxD Didx] = max(DISCR,[],2);
                    uDidx = unique(Didx);
                    for k = 1:length(uDidx)
                        locDidx = Didx == uDidx(k);
                        clusteringMatched.maxTausPerSpike(locDidx,1) = maxTaus(locDidx,uDidx(k));
                        clusteringMatched.maxTausIdxPerSpike(locDidx,1) = tauIdx(locDidx,uDidx(k));
                        clusteringMatched.maxSubTausIdxPerSpike(locDidx,1) = subTauIdx(locDidx,uDidx(k));
                    end
                    notNoiseIdx = MaxD>0;
                    clusteringMatched.ids(notNoiseIdx,1) = Didx(notNoiseIdx);
                    clusteringMatched.template_ids = unique(clusteringMatched.ids);

                    if 0 %~isempty(P.templateMatchingCut.residualPeakSmallerThanStdNoise)
                        % compute the "explained energy"
                        notNoiseIdx = find(notNoiseIdx);
                        nSmatched = length(notNoiseIdx);
                        E = zeros(size(S.spikeCut.wfs));
                        tTf = size(TupShiftDown,2);
                        for s = 1:8 %nSmatched
                            myId = clusteringMatched.ids(notNoiseIdx(s));
                            mySubTauIdx = clusteringMatched.maxSubTausIdxPerSpike(notNoiseIdx(s));
                            myTau = floor(clusteringMatched.maxTausPerSpike(notNoiseIdx(s)));
                            myStartIdx = idxstart-myTau;
                            E(notNoiseIdx(s),myStartIdx:myStartIdx+tTf-1) = squeeze(TupShiftDown(myId,:,mySubTauIdx));
                        end
                        R = S.spikeCut.wfs-E;
                        figure;
                        subplot(2,1,1);
                        plot(S.spikeCut.wfs(1:4,:)');
                        hold on
                        plot(E(1:4,:)');
                        subplot(2,1,2); plot(R(1:4,:)');
                    end

                    % shift the spikes to compute the templates
                    chk = hdsort.util.Chunker(size(S.spikeCut.wfs,1), 'chunkSize', P.spikeCutting.chunkSize, ...
                                 'chunkOverlap', 0);
                    allIDs = unique(clusteringMatched.ids);
                    nT = length(allIDs);
                    N_total_per_template = zeros(1, nT);
                    clusteringMatched.templates = zeros(nT, size(S.spikeCut.wfs,2));
                    
                    N_SPIKES_PER_TEMPLATE = 400;
                    % Select randomly a subset of spikes for each template
                    % to compute the template from
                    CUTME = [];
                    for i=1:length(allIDs)
                        myIdx = find(clusteringMatched.ids==allIDs(i), N_SPIKES_PER_TEMPLATE);
                        myID = ones(length(myIdx),1)*allIDs(i);
                        CUTME = [CUTME; myIdx(:) myID(:)];
                    end
                    wfs_temp = S.spikeCut.wfs(CUTME(:,1),:);
                    wfs_temp_aligned =hdsort.waveforms.vShift(wfs_temp, nC, -round(clusteringMatched.maxTausPerSpike(CUTME(:,1))), 1);
                        [templates_temp, N_per_template] = hdsort.util.calculateClassMeans(...
                            wfs_temp_aligned, CUTME(:,2), 'usemedian');  
                    clusteringMatched.templates  = templates_temp;
                    clusteringMatched.README = 'WARNING THE TEMPLATES AFTER MATCHING ARE NOT PROPERLY COMPUTED AT THE BOARDERS. THEY WOULD HAVE TO BE RECUT!';
                    clusteringMatched.templatesIDs = unique(clusteringMatched.ids);
                    clusteringMatched.ts = S.spikeDetectionMerged.ts(S.spikeCut.cutIdx)+clusteringMatched.maxTausPerSpike;
                    
                    %clusteringMatched.amps = S.spikeDetectionMerged.allspikes(S.spikeCut.cutIdx, 2);
                    wfstmp = hdsort.waveforms.v2t(S.spikeCut.wfs(:,:), nC);
                    [clusteringMatched.ampsChIdx, clusteringMatched.amps] = hdsort.waveforms.maxChannel(wfstmp);
                    %[clusteringMatched.amps, clusteringMatched.ampsChIdx] = max(wfstmp, [], 2);
                    clear wfstmp
                   % [amps, ampsChIdx] = max(S.spikeCut.wfs, [], 2);
                   % clusteringMatched.amps = amps;
                    
                    if 0
                        N = 1000;
                        figure; subplot(2,1,1); plot(S.spikeCut.wfs(1:N,:)');
                        subplot(2,1,2);
                        plot(clusteringMatched.spikeCutAligned(1:N,:)');
                    end
                end
            end
            saveIfDoesNotExist(S.files.botm_matching_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.botm_matching_file, 'clusteringMatched', '-v7.3'); 
            disp('Done.')             
        end
        S.clusteringMatched = clusteringMatched;
    end
    %% --------------------------------------------------------------------
    function mergeClustersMS()
        if ~P.mergeTemplates.merge
            return
        end        
        try
            load(S.files.merge_ms_clusters);
            disp('Already merged...');          
        catch   
            disp('Merging MS clusters...');
            clusteringMerged = [];
            % in S.clusteringMatched.ids the ids are going from 0 to nT with
            % 0 being not matched and noise, thus there is one "real" template less
            % than units (if any spike has an ID of 0 !!)
            realTemplateIdx = find(S.clusteringMatched.template_ids>0);
            nT = length(realTemplateIdx);
            if nT == 0
                % we have only the noise cluster
                D = [];
                maxT = [];
                groups = {};
            else
                % resample only templates that do not have id==0
                resampledTemplates = hdsort.waveforms.tResample(hdsort.waveforms.v2t(...
                    S.clusteringMatched.templates(realTemplateIdx,:), nC), 3, 1);                

                % the indices in the groups are indices into
                % realTemplateIdx which is index into S.clusteringMatched.template_ids
                [groups maxT D] = hdsort.waveforms.mergeTemplates(resampledTemplates, S.noise.meanNoiseStd,...
                    'maxRelativeDistance', P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
                    'minCorrelation', P.mergeTemplates.atCorrelation);
            end
            clusteringMerged.D = D;
            clusteringMerged.maxT = maxT;
            clusteringMerged.groups = groups;            
            clusteringMerged.templates = []; %zeros(length(groups), P.spikeCutting.Tf*nC);
            clusteringMerged.ids = zeros(length(S.clusteringMatched.ids),1);
                      
            % we need to merge now the indices in ids accoring to the groups in
            % groups. But groups points into realTemplateIdx !!
            for i=1:length(groups)
                group_ids = S.clusteringMatched.template_ids(realTemplateIdx(clusteringMerged.groups{i}));
                idx = ismember(S.clusteringMatched.ids, group_ids);
                clusteringMerged.ids(idx) = i;
                % THIS IS A HACK! IN PRINCIPLE, THE NEW TEMPLATE MUST BE A
                % WEIGHTED AVERAGE OF THE GROUP TEMPLATES !
                %% THIS DOESNT WORK AT THE MOMENT!
%                 clusteringMerged.templates(i,:) = mean(S.clusteringMatched.templates(group_ids,:),1);
            end 
            
            saveIfDoesNotExist(S.files.merge_ms_clusters);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.merge_ms_clusters, 'clusteringMerged', '-v7.3');
            fprintf('Templates after merging: %d\n', length(clusteringMerged.groups));
        end
        S.clusteringMerged = clusteringMerged;
    end    

    %% --------------------------------------------------------------------
    function runBOTM()
        if ~P.botm.run
            return
        end

        try
            load(S.files.botm_file);
            disp('already sorted with botm...');
            botm = [];
        catch
            disp('BOTM sorting...');
            T = S.clusteringMerged.templates;
            Cest = S.noise.CestS;
            Cest.CCol(1:nC,:) = Cest.CCol(1:nC,:) + diag(diag(Cest.CCol(1:nC,:)));
            Cest.CCol = Cest.CCol/2;
            Cest = hdsort.noise.Covest2(Cest);
            idxstart = 1 + P.spikeCutting.cutLeft - P.botm.cutLeft;
            idxstop  = idxstart + P.botm.Tf - 1;
            botm.templates = hdsort.waveforms.vSubsel(T, nC, idxstart:idxstop);
            botmsorter = botm.BOTM(Cest, P.botm.Tf, botm.templates,...
                        'upsample', 3, 'spikePrior', P.botm.prior,...
                        'max_num_of_channels_per_template', 30, 'adaptOnInit', 1);    
            
            botmsorter.DH = DS;    
            if 0
                plot.SliderDataAxes({DS botmsorter})
            end
            botm.gdfUnclean = botmsorter.sort(DS);
            % Remove artefacts from gdf
            botm.removeIdx = epoch.findPointsInEpochs(botm.gdfUnclean(:,2), S.artefactDetection.epochs)>0;               
            botm.gdf = botm.gdfUnclean(~botm.removeIdx,:);
            saveIfDoesNotExist(S.files.botm_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
            save(S.files.botm_file, 'botm');
            S.botmsorter = botmsorter;
            disp('Done.')    
        end
        S.botm = botm; 
    end
    %% --------------------------------------------------------------------
    function computeBOTMTemplates()
        if ~isempty(gdfbotm)
            if exist(S.files.botm_templates_file, 'file')
                disp('loading botm templates...');
                D = load(S.files.botm_templates_file);
                templatesAfterBOTM = D.templatesAfterBOTM;
                clear D
            else
                disp('computing botm templates...');
                uclasses = unique(gdfbotm(:,1));
                templatesAfterBOTM = zeros(length(uclasses), P.spike_cut_Tf*nC);
                for i=1:length(uclasses)
                    myts = gdfbotm(gdfbotm(:,1)==uclasses(i),2);
                    if size(myts,1) > 1000
                        myts = myts(1:1000);
                    end
                    spikes = DS.getWaveform(myts, P.spike_cut_cutLeft, P.spike_cut_Tf);
                    templatesAfterBOTM(i,:) = mean(spikes,1);
                end
                saveIfDoesNotExist(S.files.botm_templates_file);                 
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
                save(S.files.botm_templates_file, 'templatesAfterBOTM');
                disp('Done.')    
            end            
        end
        S.botmDetails = botmDetails; clear botmDetails;
    end
    %% --------------------------------------------------------------------
    function residualArtefactDetection()
        if ~isempty(gdfbotm)
            if exist(S.files.resArtefact_file, 'file')
                disp('residual artefacts already detected');
                D = load(S.files.resArtefact_file);
                gdfbotmcleaned = D.gdfbotmcleaned;
                resArtefactEpochs = D.resArtefactEpochs; clear D;
            else
                disp('Residual atrefact detection ...');
                SpSoBOTM = spiketrain.SpikeSortingContainer('botm', gdfbotm, ...
                    'templateCutLeft', S.P.spike_cut_cutLeft, ...
                    'templateCutLength', S.P.spike_cut_Tf, ...
                    'templateWfs', hdsort.waveforms.v2t(templatesAfterBOTM, nC),...
                    'wfDataSource', DS, 'nMaxSpikesForTemplateCalc', 1000);
                DS.addSpikeSorting(SpSoBOTM); DS.setActiveSpikeSortingIdx('botm')            
                DS.bReturnSortingResiduals = 1;
                A = ana.moritzheimdahl.ResidualArtefactDetector(DS);

                resArtefactEpochs = A.getEpochs();
                gdfbotmcleaned = gdfbotm;
                removeIdx = epoch.findPointsInEpochs(gdfbotm(:,2), resArtefactEpochs)>0;
                gdfbotmcleaned(removeIdx,:) = [];      
                saveIfDoesNotExist(S.files.resArtefact_file);                 
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end                 
                save(S.files.resArtefact_file, 'resArtefactEpochs', 'gdfbotmcleaned');
                DS.bReturnSortingResiduals = 0;
                disp('Done.')    
            end
        end
        S.botmResiduals = botmResiduals; clear botmResiduals;
    end
    %% --------------------------------------------------------------------
    function computeBOTMCleanedTemplates()
        if ~isempty(gdfbotmcleaned)
            if exist(S.files.botm_templates_cleaned_file, 'file')
                disp('loading botm cleaned templates...');
                D = load(S.files.botm_templates_cleaned_file);
                templatesAfterBOTMcleaned = D.templatesAfterBOTMcleaned;
                clear D
            else
                disp('computing botm cleaned templates...');
                uclasses = unique(gdfbotmcleaned(:,1));
                templatesAfterBOTMcleaned = zeros(length(uclasses), P.spike_cut_Tf*nC);
                for i=1:length(uclasses)
                    myts = gdfbotmcleaned(gdfbotmcleaned(:,1)==uclasses(i),2);
                    if size(myts,1) > 1000
                        myts = myts(1:1000);
                    end
                    spikes = DS.getWaveform(myts, P.spike_cut_cutLeft, P.spike_cut_Tf);
                    templatesAfterBOTMcleaned(i,:) = mean(spikes,1);
                end
                saveIfDoesNotExist(S.files.botm_templates_cleaned_file);
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                save(S.files.botm_templates_cleaned_file, 'templatesAfterBOTMcleaned');
                disp('Done.')    
            end            
        end
        S.botmCleaned = botmCleaned; clear botmCleaned;                
    end

    %% --------------------------------------------------------------------
    function saveIfDoesNotExist(fname)
        % This function is only needed since we have the sorter
        % parallelized and it can happen that the very same sorting runs
        % twice on the same files. In this case, the first job to finish
        % should be the one that keeps the files, they should not be
        % overwritten by the slower one.
        if exist(fname, 'file')>0
            warning(sprintf('File: %s does already exist, WILL NOT BE OVERWRITTEN - sorting ended!', fname));
            S.STOP_ME_BECAUSE_I_AM_SLOW = true;
            %S.S.STOP_ME_BECAUSE_I_AM_SLOW = S.STOP_ME_BECAUSE_I_AM_SLOW;
            % Wait a bit because of the file system issue...
            pause(10); 
            return
        end
    end
%     %% --------------------------------------------------------------------
%     function mergeClustersBotm()
%         if 1 %isempty(gdfbotmcleaned)
%             return
%         end
%     end
end