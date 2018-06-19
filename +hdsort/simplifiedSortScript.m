function [S P] = simplifiedSortScript(WFS, noise, dpath, name, varargin)
    S.STOP_ME_BECAUSE_I_AM_SLOW = false;
    
    P.spikeDetection = [];
    P.artefactDetection = [];
    P.botm = [];
    P.spikeCutting = [];
    
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
    P.clustering.minSpikesPerCluster = 10;
    
    % Template Matching on Cut Spikes
    P.templateMatchingCut.chunkSize = 10000;
    
    P.templateMatchingCut.prior = .0001;
    P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];
    
    % Merge Templates after template matching
    P.mergeTemplates.merge = 1;
	P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .96;
    P.mergeTemplates.ifMaxDistSmaller = 2.5;      % in units std of noise
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 25;

    % Used later for full template estimation:
    P.templateEstimation.cutLeft = 10;
    P.templateEstimation.Tf = 55;
    P.templateEstimation.maxSpikes = 600;
    
    P = hdsort.util.parseInputs(P, varargin, 'error');
    
    % Checks
    assert(P.featureExtraction.Tf < P.spikeAlignment.Tf, 'Feature extraction spikeforms must be shorter than the alignment spikeforms!')
    assert(P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft >= 0, 'Connot prewhiten spikes with this cutleft!');
    % Save header
    readme = 'The files in this folder were created by the simplifiedSortScript.m script.';
    
    % Define global names
    S.P = P;
    S.name = name;
    S.dpath = dpath;
    S.dprefix = fullfile(dpath, name);
    
    if ~exist(S.dpath, 'file')
        fprintf('Output directory does not exists. Trying to create...')
        mkdir(S.dpath);
        fprintf(' success.\n')
    end   
    
    % Define file names
    S.files.spike_cut_file         = [S.dprefix '.040spikes_cut.mat'];
    S.files.spike_aligned_file     = [S.dprefix '.050spikes_aligned.mat'];
    S.files.cov_file               = [S.dprefix '.060cov.mat'];
    S.files.prewh_spike_file       = [S.dprefix '.070spikes_prewhitened.mat'];
    S.files.fet_spike_file         = [S.dprefix '.080spikes_features.mat'];
    S.files.meanshift_spike_file   = [S.dprefix '.090clusters_meanshift.mat'];
    S.files.botm_matching_file     = [S.dprefix '.100botm_matching.mat'];    
    S.files.merge_ms_clusters      = [S.dprefix '.110clusters_meanshift_merged.mat'];     
    save([fullfile(dpath, name) '.P.mat'], 'S', 'readme');
    
    % RUN
    prepareInput();         if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    alignSpikes();          if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    prewhitenSpikes();      if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    featureExtraction();    if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runMeanShift();         if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    runBOTMMatching();      if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    mergeClustersMS();      if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
    
    disp('############################');
    disp(['Sorting of ' S.name ' finished.']);
    
    %----------------------------------------------------------------------
    % Detected sp:  +++++ ++ +++++++++++++ +++++++ ++++++++++ ++++++ +++++
    % Cut spikes :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    % Aligned sp :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Prewhitend :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Features   :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Clustered  :     +     +           +     + +                 + +
    % Matched    :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    
    %% --------------------------------------------------------------------
    function prepareInput()
        S.srate = WFS.samplesPerSecond;
        S.nC = WFS.getNChannels();
        S.spikeCut.wfs = WFS;
        nSP = WFS.getNSamples();
        
        if P.spikeAlignment.maxSpikes < nSP
            alignIdx = randperm(nSP);
            S.spikeCut.alignIdx = sort(alignIdx(1:P.spikeAlignment.maxSpikes));
        else
            S.spikeCut.alignIdx = 1:nSP;
        end
        S.spikeCut.nSpikesToAlign = length(S.spikeCut.alignIdx);
        S.spikeCut.unalignedwfs = S.spikeCut.wfs(S.spikeCut.alignIdx, :);
        S.spikeCut.cutIdx = 1:nSP;
        
        % Override this variable (because this parameter was set before):
        P.spikeCutting.Tf = size(S.spikeCut.wfs, 2) / S.nC;
        P.spikeCutting.cutLeft = WFS.getCutLeft();
        
        % S.noise:
        S.noise = noise;
        maxTf_ = max([P.featureExtraction.Tf P.spikeAlignment.Tf]);
        S.noise.C_time_aligned = hdsort.noise.ccol2Cte(S.noise.CestS.CCol, maxTf_);
        
        % S.spikeDetectionMerged:
        if nSP
            gdf = WFS.getGdf();
            S.spikeDetectionMerged.ts = gdf(:,2);
            clear gdf
        else
            S.spikeDetectionMerged.ts = [];
        end
        assert( all(S.spikeDetectionMerged.ts > 0), '!!!')
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
                t_align = tic;
                
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
                        S.nC, idxstart:idxstop));    

                    spikeAligned.tau = zeros(size(spikeAligned.unalignedwfs,1),1);
                    spikeAligned.maxIdx = P.spikeAlignment.cutLeft;
                    spikeAligned.idxstart = idxstart;
                    spikeAligned.idxstop = idxstop;

                    spikeAligned.restrictToIdx = spikeAligned.maxIdx-2:spikeAligned.maxIdx+5;
                    tic
                    if strcmp(P.spikeAlignment.method, 'onMax')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnMax(spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMean')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnUpsampleMean(spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', P.spikeAlignment.maxIdx,...
                            'maxShiftPerIter', 3,...
                            'maxIter', P.spikeAlignment.maxIterations, ...
                            'initAlignment', P.spikeAlignment.initAlignment);   
%                             'restrictToIdx',  spikeAligned.restrictToIdx,...                        
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.util.alignWaveformsUpsampleMax(spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                    elseif strcmp(P.spikeAlignment.method, 'onUpsampledMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.alignWaveformsUpsampleMax(-spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2); 
                        spikeAligned.wfs = -spikeAligned.unalignedwfs;
                    elseif strcmp(P.spikeAlignment.method, 'onMin')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            hdsort.waveforms.vAlignOnMax(-spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(P.spikeAlignment.method, 'onAverageMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.vAlignOnAverageMaxSample(spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(P.spikeAlignment.method, 'onAverageMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            hdsort.waveforms.vAlignOnAverageMaxSample(-spikeAligned.unalignedwfs, S.nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(P.spikeAlignment.method, 'none')
                        spikeAligned.wfs = spikeAligned.unalignedwfs;
                    else
                        error(['unkown alignement method' P.spikeAlignment.method]);
                    end
                    toc
                end
                % cut away alignement artefacts at the ends
                saveIfDoesNotExist(S.files.spike_aligned_file);                 
                if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                
                spikeAligned.time = toc(t_align);
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
            t_prewhiten = tic;
            % Cut the aligned waveforms to get rid of alignement artefacts
            % and reduce to the final waveform that will be used for
            % feature extraction
                idxstart = 1 + P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft;
                idxstop  = idxstart + P.featureExtraction.Tf - 1;   
                ccol_loaded = S.noise.CestS.CCol;
                ccol_loaded(1:S.nC, 1:S.nC) = ccol_loaded(1:S.nC, 1:S.nC) + eye(S.nC) * .1 * mean(diag(ccol_loaded(1:S.nC, 1:S.nC)));                
                spikePrewhitened.ccol_loaded = ccol_loaded/2;
                spikePrewhitened.C = hdsort.noise.ccol2Cte(spikePrewhitened.ccol_loaded, P.featureExtraction.Tf);
                spikePrewhitened.U = chol(spikePrewhitened.C);            
            if isempty(S.spikeAligned.wfs)
                % Prewhiten
                spikePrewhitened.wfs = [];
            else
                spikePrewhitened.wfs = hdsort.waveforms.vSubsel(S.spikeAligned.wfs, S.nC, idxstart:idxstop);      
                % Build the noise covariance matrix and load it
                spikePrewhitened.wfs = spikePrewhitened.wfs/spikePrewhitened.U; 
            end
            saveIfDoesNotExist(S.files.prewh_spike_file);                 
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            spikePrewhitened.time = toc(t_prewhiten);
            save(S.files.prewh_spike_file, 'spikePrewhitened', '-v7.3');
            disp('Done.')    
        end
        S.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
    end

    %% --------------------------------------------------------------------
    function featureExtraction()
        try
            spikeFeatures = [];
            load(S.files.fet_spike_file );  
            disp('Features already calculated...');
        catch
            disp('Calculating features...');
            t_fet = tic;
            if ~isempty(S.spikePrewhitened.wfs)
                [spikeFeatures.X, spikeFeatures.PCs] = hdsort.util.dimReductionPCA(S.spikePrewhitened.wfs,...
                    P.featureExtraction.nDims, [], 3*1000000);
            else
                spikeFeatures.X = [];
            end
            saveIfDoesNotExist(S.files.fet_spike_file);
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            spikeFeatures.time = toc(t_fet);
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
            t_meanShift = tic;
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
            clustering.time = toc(t_meanShift);
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
            t_BOTMMatching = tic;
            nSpCut = size(S.spikeCut.wfs(:,:),1);
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
                    C = (S.noise.C_time_aligned + diag(diag(S.noise.C_time_aligned)))/2;

                    % Do matching only on temporal window defined by the aligned
                    % templates
                    idxstart = 1 + P.spikeCutting.cutLeft - P.spikeAlignment.cutLeft;
                    
                    [D maxTaus TupShiftDown FupShiftDown tauRange tauIdx subTauRange subTauIdx] = hdsort.waveforms.vTemplateMatching(...
                        S.spikeCut.wfs, T, S.nC, idxstart, 'maxShift', 5, 'upsample', 5, 'noiseCovariance', C, 'chunkSize', P.templateMatchingCut.chunkSize);
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
                    
                    % shift the spikes to compute the templates
                    chk = hdsort.util.Chunker(size(S.spikeCut.wfs,1), 'chunkSize', P.templateMatchingCut.chunkSize, ...
                                 'chunkOverlap', 0);
                    allIDs = unique(clusteringMatched.ids);
                    nT = length(allIDs);
                    N_total_per_template = zeros(1, nT);
                    clusteringMatched.templates = zeros(nT, size(S.spikeCut.wfs,2));
                    
                    N_SPIKES_PER_TEMPLATE = 400;
                    % Select randomly a subset of spikes for each template
                    % to compute the template from
                    nSpikesPerTemplate = zeros(length(allIDs),1);
                    CUTME = [];
                    for i=1:length(allIDs)
                        myIdx = find(clusteringMatched.ids==allIDs(i), N_SPIKES_PER_TEMPLATE);
                        myID = ones(length(myIdx),1)*allIDs(i);
                        CUTME = [CUTME; myIdx(:) myID(:)];
                        nSpikesPerTemplate(i) = length(myIdx);
                    end
                    wfs_temp = S.spikeCut.wfs(CUTME(:,1),:);
                    wfs_temp_aligned = hdsort.waveforms.vShift(wfs_temp, S.nC, -round(clusteringMatched.maxTausPerSpike(CUTME(:,1))), 1);
                        [templates_temp, N_per_template] = hdsort.util.calculateClassMeans(...
                            wfs_temp_aligned, CUTME(:,2), 'usemedian');  
                    clusteringMatched.templates  = templates_temp;
                    clusteringMatched.nSpikesPerTemplate = nSpikesPerTemplate;
                    clusteringMatched.README = 'WARNING THE TEMPLATES AFTER MATCHING ARE NOT PROPERLY COMPUTED AT THE BORDERS. THEY SHOULD BE RECUT AT THE EDGES!';
                    clusteringMatched.templatesIDs = unique(clusteringMatched.ids);
                    clusteringMatched.ts = S.spikeDetectionMerged.ts(S.spikeCut.cutIdx)+clusteringMatched.maxTausPerSpike;
                end
                
                wfstmp = hdsort.waveforms.v2t(S.spikeCut.wfs(:,:), S.nC);
                [clusteringMatched.ampsChIdx, clusteringMatched.amps] = hdsort.waveforms.maxChannel(wfstmp);
                clear wfstmp
            end
            saveIfDoesNotExist(S.files.botm_matching_file);
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            clusteringMatched.time = toc(t_BOTMMatching);
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
            t_merge = tic;
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
                    S.clusteringMatched.templates(realTemplateIdx,:), S.nC), 3, 1);                

                % the indices in the groups are indices into
                % realTemplateIdx which is index into S.clusteringMatched.template_ids
                [groups maxT D] = hdsort.waveforms.mergeTemplates(resampledTemplates, S.noise.meanNoiseStd,...
                    'maxRelativeDistance', P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
                    'minCorrelation', P.mergeTemplates.atCorrelation);
                
            end
            clusteringMerged.D = D;
            clusteringMerged.maxT = maxT;
            clusteringMerged.groups = groups;            
            clusteringMerged.templates = zeros(length(groups), P.spikeCutting.Tf*S.nC);
            clusteringMerged.ids = zeros(length(S.clusteringMatched.ids),1);
                      
            % now, we need to merge the indices in ids accoring to the groups in
            % groups. But groups points into realTemplateIdx !!
            for i=1:length(groups)
                group_idx = realTemplateIdx(clusteringMerged.groups{i});
                group_ids = S.clusteringMatched.template_ids(group_idx);
                idx = ismember(S.clusteringMatched.ids, group_ids);
                % get number of spikes for each template for weighted
                % average
                nSpikesPerTemplate = zeros(1,length(group_ids));
                for kk=1:length(group_ids)
                    nSpikesPerTemplate(kk) = sum(S.clusteringMatched.ids==group_ids(kk));
                end
                clusteringMerged.ids(idx) = i;
                clusteringMerged.templates(i,:) = nSpikesPerTemplate*S.clusteringMatched.templates(group_idx,:)/sum(nSpikesPerTemplate);
                if 0
                    %%
                    figure
                    subplot(2,1,1)
                    plot(S.clusteringMatched.templates(group_idx,:)');
                    subplot(2,1,2)
                    plot(clusteringMerged.templates(i,:));
                end
            end 
            
            saveIfDoesNotExist(S.files.merge_ms_clusters);
            if S.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            clusteringMerged.time = toc(t_merge);
            save(S.files.merge_ms_clusters, 'clusteringMerged', '-v7.3');
            fprintf('Templates after merging: %d\n', length(clusteringMerged.groups));
        end
        S.clusteringMerged = clusteringMerged;
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


end