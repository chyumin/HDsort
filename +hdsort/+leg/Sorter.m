classdef Sorter < handle
    
    properties
        P
        path
        name
        prefix
        files
        LEGindices % is a list with the channels in this LEG
        
        buffer
        WFS
        noise
        
        STOP_ME_BECAUSE_I_AM_SLOW
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Sorter(WFS, noise, sortingPath, sortingName, LEGindices, varargin)
            self.path = sortingPath;
            self.name = sortingName;
            self.prefix = fullfile(self.path, self.name);
            self.LEGindices = LEGindices;
            
            self.STOP_ME_BECAUSE_I_AM_SLOW = false;
            self.WFS = WFS;
            self.noise = noise;
            
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
            P.templateMatchingCut.prior = .0001;
            P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];
            P.templateMatchingCut.chunkSize = 5000;
            
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
            
            P = hdsort.util.parseInputs(P, varargin, 'split');
            
            % Checks
            assert(self.WFS.getCutLeft() >= P.spikeAlignment.cutLeft, 'Cannot align further left than cut !');
            assert(P.featureExtraction.Tf < P.spikeAlignment.Tf, 'Feature extraction spikeforms must be shorter than the alignment spikeforms!')
            assert(P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft >= 0, 'Connot prewhiten spikes with this cutLeft!');
            % Save header
            readme = 'The files in this folder were created by hdsort.leg.Sorter.';
            
            % Create sorting path:
            if ~exist(self.path, 'file')
                fprintf('Output directory does not exists. Trying to create...')
                mkdir(self.path);
                fprintf(' success.\n')
            end
            
            % Define file names
            self.files.P                      = [self.prefix '.P.mat'];
            self.files.artefact_file          = [self.prefix '.010artefacts.mat'];
            self.files.spike_det_file         = [self.prefix '.020spikes_det.mat'];
            self.files.spike_det_merged_file  = [self.prefix '.030spikes_det_merged.mat'];
            self.files.spike_cut_file         = [self.prefix '.040spikes_cut.mat'];
            self.files.spike_aligned_file     = [self.prefix '.050spikes_aligned.mat'];
            self.files.cov_file               = [self.prefix '.060cov.mat'];
            self.files.prewh_spike_file       = [self.prefix '.070spikes_prewhitened.mat'];
            self.files.fet_spike_file         = [self.prefix '.080spikes_features.mat'];
            self.files.meanshift_spike_file   = [self.prefix '.090clusters_meanshift.mat'];
            self.files.botm_matching_file     = [self.prefix '.100botm_matching.mat'];
            self.files.merge_clusters         = [self.prefix '.110clusters_meanshift_merged.mat'];
            self.files.gdf                    = [self.prefix '.gdf.mat'];
            self.files.full_template          = [self.prefix '_templates.mat'];
            
            save(self.files.P, 'P', 'readme');
            self.P = P;
            clear P
            
            % Prepare cutSpikes:
            self.buffer.srate = self.WFS.samplesPerSecond;
            self.buffer.nC = self.WFS.getNChannels();
            assert( numel(self.LEGindices) == self.buffer.nC, 'Number of channels must be the same!')
            
            nSP = self.WFS.getNSamples();
            if self.P.spikeAlignment.maxSpikes < nSP
                alignIdx = randperm(nSP);
                self.buffer.spikeCut.alignIdx = sort(alignIdx(1:self.P.spikeAlignment.maxSpikes));
            else
                self.buffer.spikeCut.alignIdx = 1:nSP;
            end
            self.buffer.spikeCut.nSpikesToAlign = length(self.buffer.spikeCut.alignIdx);
            
            disp('Load waveforms to memory...')
            self.buffer.spikeCut.unalignedwfs = self.WFS(self.buffer.spikeCut.alignIdx, :);
            disp('Loading waveforms done.')
            
            self.buffer.spikeCut.cutIdx = 1:nSP;
            
            % Prepare  spikeDetectionMerged
            if nSP
                gdf = self.WFS.getGdf();
                self.buffer.spikeDetectionMerged.ts = gdf(:,2);
                clear gdf
            else
                self.buffer.spikeDetectionMerged.ts = [];
            end
            assert( all(self.buffer.spikeDetectionMerged.ts > 0), '!!!')
            
            % Prepare noise
            maxTf_ = max([self.P.featureExtraction.Tf self.P.spikeAlignment.Tf]);
            self.noise.C_time_aligned = hdsort.noise.ccol2Cte(self.noise.CestS.CCol, maxTf_);
            
        end
        
        %% --------------------------------------------------------------------
        function S = runAll(self, DS)
            % RUN
            self.alignSpikes();       if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.prewhitenSpikes();   if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.featureExtraction(); if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.meanShift();         if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.BOTMMatching();      if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.mergeClusters();     if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            self.saveGDF();           if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            
            if nargin > 1
                
                % Try this to reduce memory usage:
                gdf = self.buffer.gdf;
                self.buffer = [];
                self.buffer.gdf = gdf;
                clear gdf;
                
                self.estimateFullTemplate(DS);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
            end
            
            disp('############################');
            disp(['Sorting of ' self.name ' finished.']);
            
            S = self.buffer;
            S.P = self.P;
        end
        
        %% --------------------------------------------------------------------
        function alignSpikes(self)
            
            % Make sure to align the spike maximall in the window of the dead
            % time of the spike detector. Otherwise, we detect to spikes which
            % are close to each other and then we shift them to effectively
            % detect the same spike twice.
            if ~isempty(self.P.spikeAlignment.method);
                try
                    spikeAligned = [];
                    load(self.files.spike_aligned_file);
                    disp('Spikes already aligned');
                catch
                    disp('Aligning spikes...');
                    t_align = tic;
                    
                    idxstart = 1 + self.WFS.getCutLeft() - self.P.spikeAlignment.cutLeft;
                    idxstop  = idxstart + self.P.spikeAlignment.Tf - 1;
                    nSpCut = size(self.WFS,1);
                    
                    %% %%%%%%%%%%%%%%%%%
                    if nSpCut == 0
                        spikeAligned.tau = [];
                        spikeAligned.maxIdx = [];
                        spikeAligned.restrictToIdx =[];
                        spikeAligned.wfs = [];
                        spikeAligned.alignIdx = [];
                    else
                        spikeAligned.alignIdx = self.buffer.spikeCut.alignIdx;
                        spikeAligned.unalignedwfs = double(hdsort.waveforms.vSubsel(self.buffer.spikeCut.unalignedwfs,...
                            self.buffer.nC, idxstart:idxstop));
                        
                        % Free memory:
                        self.buffer.spikeCut.unalignedwfs = [];
                        
                        spikeAligned.tau = zeros(size(spikeAligned.unalignedwfs,1),1);
                        spikeAligned.maxIdx = self.P.spikeAlignment.cutLeft;
                        spikeAligned.idxstart = idxstart;
                        spikeAligned.idxstop = idxstop;
                        
                        spikeAligned.restrictToIdx = spikeAligned.maxIdx-2:spikeAligned.maxIdx+5;
                        tic
                        if strcmp(self.P.spikeAlignment.method, 'onMax')
                            [spikeAligned.wfs spikeAligned.tau] = ...
                                hdsort.waveforms.vAlignOnMax(spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMean')
                            [spikeAligned.wfs spikeAligned.tau] = ...
                                hdsort.waveforms.vAlignOnUpsampleMean(spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', self.P.spikeAlignment.maxIdx,...
                                'maxShiftPerIter', 3,...
                                'maxIter', self.P.spikeAlignment.maxIterations, ...
                                'initAlignment', self.P.spikeAlignment.initAlignment);
                            %                             'restrictToIdx',  spikeAligned.restrictToIdx,...
                        elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMax')
                            [spikeAligned.tau spikeAligned.wfs] = ...
                                hdsort.util.alignWaveformsUpsampleMax(spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                        elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMin')
                            [spikeAligned.tau spikeAligned.wfs] = ...
                                hdsort.waveforms.alignWaveformsUpsampleMax(-spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                            spikeAligned.wfs = -spikeAligned.unalignedwfs;
                        elseif strcmp(self.P.spikeAlignment.method, 'onMin')
                            [spikeAligned.wfs spikeAligned.tau] = ...
                                hdsort.waveforms.vAlignOnMax(-spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                            spikeAligned.wfs = -spikeAligned.wfs;
                        elseif strcmp(self.P.spikeAlignment.method, 'onAverageMax')
                            [spikeAligned.tau spikeAligned.wfs] = ...
                                hdsort.waveforms.vAlignOnAverageMaxSample(spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        elseif strcmp(self.P.spikeAlignment.method, 'onAverageMin')
                            [spikeAligned.tau spikeAligned.wfs] = ...
                                hdsort.waveforms.vAlignOnAverageMaxSample(-spikeAligned.unalignedwfs, self.buffer.nC,...
                                'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                            spikeAligned.wfs = -spikeAligned.wfs;
                        elseif strcmp(self.P.spikeAlignment.method, 'none')
                            spikeAligned.wfs = spikeAligned.unalignedwfs;
                        else
                            error(['unkown alignement method' self.P.spikeAlignment.method]);
                        end
                        toc
                    end
                    
                    self.saveIfDoesNotExist(self.files.spike_aligned_file);
                    if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                    
                    spikeAligned.time = toc(t_align);
                    save(self.files.spike_aligned_file, 'spikeAligned', '-v7.3');
                    disp('Aligning spikes done.')
                end
            end
            self.buffer.spikeAligned = spikeAligned; clear spikeAligned;
            self.buffer.spikeAligned.unalignedwfs = []; % free memory
            self.buffer.spikeCut.unalignedwfs = []; % free memory
        end
        
        %% --------------------------------------------------------------------
        function prewhitenSpikes(self)
            try
                spikePrewhitened = [];
                load(self.files.prewh_spike_file);
                disp('Spikes already prewhitened...');
            catch
                disp('Prewhitening spikes...');
                t_prewhiten = tic;
                % Cut the aligned waveforms to get rid of alignement artefacts
                % and reduce to the final waveform that will be used for
                % feature extraction
                idxstart = 1 + self.P.spikeAlignment.cutLeft - self.P.featureExtraction.cutLeft;
                idxstop  = idxstart + self.P.featureExtraction.Tf - 1;
                ccol_loaded = self.noise.CestS.CCol;
                ccol_loaded(1:self.buffer.nC, 1:self.buffer.nC) = ccol_loaded(1:self.buffer.nC, 1:self.buffer.nC) + eye(self.buffer.nC) * .1 * mean(diag(ccol_loaded(1:self.buffer.nC, 1:self.buffer.nC)));
                spikePrewhitened.ccol_loaded = ccol_loaded/2;
                spikePrewhitened.C = hdsort.noise.ccol2Cte(spikePrewhitened.ccol_loaded, self.P.featureExtraction.Tf);
                spikePrewhitened.U = chol(spikePrewhitened.C);
                if isempty(self.buffer.spikeAligned.wfs)
                    % Prewhiten
                    spikePrewhitened.wfs = [];
                else
                    spikePrewhitened.wfs = hdsort.waveforms.vSubsel(self.buffer.spikeAligned.wfs, self.buffer.nC, idxstart:idxstop);
                    % Build the noise covariance matrix and load it
                    spikePrewhitened.wfs = spikePrewhitened.wfs/spikePrewhitened.U;
                end
                self.saveIfDoesNotExist(self.files.prewh_spike_file);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                spikePrewhitened.time = toc(t_prewhiten);
                save(self.files.prewh_spike_file, 'spikePrewhitened', '-v7.3');
                disp('Prewhitening spikes done.')
            end
            self.buffer.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
        end
        
        %% --------------------------------------------------------------------
        function featureExtraction(self)
            try
                spikeFeatures = [];
                load(self.files.fet_spike_file );
                disp('Features already calculated...');
            catch
                disp('Calculating features...');
                t_fet = tic;
                if ~isempty(self.buffer.spikePrewhitened.wfs)
                    [spikeFeatures.X, spikeFeatures.PCs] = hdsort.util.dimReductionPCA(self.buffer.spikePrewhitened.wfs,...
                        self.P.featureExtraction.nDims, [], 3*1000000);
                else
                    spikeFeatures.X = [];
                end
                self.saveIfDoesNotExist(self.files.fet_spike_file);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                spikeFeatures.time = toc(t_fet);
                save(self.files.fet_spike_file, 'spikeFeatures', '-v7.3');
                disp('Prewhitening spikes done.')
            end
            self.buffer.spikeFeatures = spikeFeatures; clear spikeFeatures;
            self.buffer.spikePrewhitened.wfs = []; % for memory minimization
        end
        
        %% --------------------------------------------------------------------
        function meanShift(self)
            try
                clustering = [];
                load(self.files.meanshift_spike_file);
                disp('Already clustered...');
            catch
                disp('Mean Shift Clustering...');
                t_meanShift = tic;
                if ~isempty(self.P.clustering.meanShiftBandWidth)
                    clustering.bandwidth = self.P.clustering.meanShiftBandWidth;
                else
                    clustering.bandwidth = sqrt(self.P.clustering.meanShiftBandWidthFactor*self.P.featureExtraction.nDims);
                end
                nSpFet = size(self.buffer.spikeFeatures.X,1);
                self.buffer.maxNBandwidthIncreases = 1;
                self.buffer.bandwidthIncreaseFactor = 1.3;
                if ~isempty(self.buffer.spikeFeatures)
                    % Select a random subset of spikes:
                    if self.P.clustering.maxSpikes < nSpFet
                        clusterIdx = randperm(nSpFet);
                        clustering.clusterIdx = sort(clusterIdx(1:self.P.clustering.maxSpikes));
                    else
                        clustering.clusterIdx = 1:nSpFet;
                    end
                    clusteringFeatures = self.buffer.spikeFeatures.X(clustering.clusterIdx,:)';
                    
                    tBeforeMeanShift = tic;
                    [~, ~, clustMembsCell] = meanshift.MeanShiftClusterIncreaseBW(...
                        clusteringFeatures, clustering.bandwidth, 0, self.P.clustering.minSpikesPerCluster, ...
                        self.buffer.maxNBandwidthIncreases, self.buffer.bandwidthIncreaseFactor);
                    tMeanShift = toc(tBeforeMeanShift);
                    fprintf('Mean Shift Clustering took %.1f sec\n', tMeanShift);
                    
                    [clustering.ids, clustering.clusterCenter] = meanshift.MeanShiftClusterBundleResult(...
                        clusteringFeatures', clustMembsCell, self.P.clustering.minSpikesPerCluster);
                    
                    %*** clustering.ids==0 --> noise cluster ***
                    clustering.classes = unique(clustering.ids);
                    
                    % Compute the mean waveforms:
                    clusteredAlignedSpikes = self.buffer.spikeAligned.wfs(clustering.clusterIdx,:);
                    clustering.templatesAligned = hdsort.util.calculateClassMeans(...
                        clusteredAlignedSpikes, clustering.ids, 'usemedian');
                    clear clusteredAlignedSpikes
                    
                    clusteredCutSpikes = self.WFS(self.buffer.spikeAligned.alignIdx(clustering.clusterIdx),:);                    
                    clustering.templatesCut = hdsort.util.calculateClassMeans(...
                        clusteredCutSpikes, clustering.ids, 'usemedian'); % templates in order of unique(clustering.ids)
                    clear clusteredCutSpikes
                    
                else
                    clustering.ids = [];
                    clustering.clusterCenter = [];
                    clustering.classes = [];
                    clustering.templatesCut  = [];
                    clustering.templatesAligned  = [];
                end
                self.saveIfDoesNotExist(self.files.meanshift_spike_file);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                clustering.time = toc(t_meanShift);
                save(self.files.meanshift_spike_file, 'clustering', '-v7.3');
            end
            self.buffer.clustering = clustering; clear clustering;
            self.buffer.spikeAligned.wfs = []; % for memory minimization
        end
        
        %% --------------------------------------------------------------------
        function BOTMMatching(self)
            
            recalc = false;
            try
                clusteringMatched = [];
                load(self.files.botm_matching_file);
                disp('BOTM matched with lda...');
            catch
                recalc = true;
            end
            
            if recalc
                disp('BOTM matching...');
                t_BOTMMatching = tic;
                nSpCut = size(self.WFS(:,:),1);
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
                    clusteringMatched.template_ids = [];
                    clusteringMatched.templatesIDs  = [];
                    clusteringMatched.amplitudeChannelIdx = [];
                    clusteringMatched.amplitudes = [];
                else
                    clusteringMatched.ids = zeros(nSpCut,1);
                    clusteringMatched.template_ids = 0;
                    clusteringMatched.spikeCutAligned = [];
                    clusteringMatched.templates = [];
                    clusteringMatched.maxTausPerSpike = zeros(nSpCut,1);
                    clusteringMatched.ts = self.buffer.spikeDetectionMerged.ts(self.buffer.spikeCut.cutIdx);
                    clusteringMatched.amplitudeChannelIdx = ones(nSpCut,1);
                    clusteringMatched.amplitudes = zeros(nSpCut,1);
                    
                    if size(self.buffer.clustering.templatesCut,1) < 2
                        disp('Cannot run BOTM Matching, not enough templates!');
                    else
                        T = self.buffer.clustering.templatesAligned(2:end,:); % remove the noise template
                        C = (self.noise.C_time_aligned + diag(diag(self.noise.C_time_aligned)))/2;
                        
                        % Do matching only on temporal window defined by the aligned
                        % templates
                        spikeCuttingCutLeft = self.WFS.getCutLeft();
                        idxstart = 1 + spikeCuttingCutLeft - self.P.spikeAlignment.cutLeft;
                        
                        [D, maxTaus, TupShiftDown, FupShiftDown, tauRange, tauIdx, ...
                            subTauRange, subTauIdx, amplitudeChannelIdx, amplitudes] = ...
                            hdsort.waveforms.vTemplateMatching(self.WFS, T, self.buffer.nC, idxstart, ...
                            'maxShift', 5, 'upsample', 5, 'noiseCovariance', C, ...
                            'chunkSize', self.P.templateMatchingCut.chunkSize);
                        
                        assert(any(D(:)~=0.0), 'SpikeCut file corrupt!')
                        
                        F = FupShiftDown(:,:,1);
                        EF = diag(T*F');                     % compute energies
                        Prior  = self.P.templateMatchingCut.prior;
                        DISCR = D - .5 * repmat(EF', nSpCut, 1)  + log(Prior);  % compute botm dicriminant
                        clusteringMatched.maxTausPerSpikeAndFilter = maxTaus;
                        clusteringMatched.tauRange = tauRange;
                        clusteringMatched.subTauRange = subTauRange;
                        clusteringMatched.tauIdx = tauIdx;
                        clusteringMatched.subTauIdx = subTauIdx;
                        clusteringMatched.TupShiftDown = TupShiftDown;
                        clusteringMatched.FupShiftDown = FupShiftDown;
                        clusteringMatched.amplitudeChannelIdx = amplitudeChannelIdx;
                        clusteringMatched.amplitudes = amplitudes;
                        
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
                        allIDs = unique(clusteringMatched.ids);
                        nT = length(allIDs);
                        N_total_per_template = zeros(1, nT);
                        clusteringMatched.templates = zeros(nT, size(self.WFS,2));
                        
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
                        wfs_temp = self.WFS(CUTME(:,1),:);
                        wfs_temp_aligned = hdsort.waveforms.vShift(wfs_temp, self.buffer.nC, -round(clusteringMatched.maxTausPerSpike(CUTME(:,1))), 1);
                        [templates_temp, N_per_template] = hdsort.util.calculateClassMeans(...
                            wfs_temp_aligned, CUTME(:,2), 'usemedian');
                        clusteringMatched.templates  = templates_temp;
                        clusteringMatched.nSpikesPerTemplate = nSpikesPerTemplate;
                        
                        clusteringMatched.README = 'WARNING THE TEMPLATES AFTER MATCHING ARE NOT PROPERLY COMPUTED AT THE BORDERself.buffer. THEY SHOULD BE RECUT AT THE EDGES!';
                        clusteringMatched.templatesIDs = unique(clusteringMatched.ids);
                        clusteringMatched.ts = self.buffer.spikeDetectionMerged.ts(self.buffer.spikeCut.cutIdx)+clusteringMatched.maxTausPerSpike;
                    end
                end
                self.saveIfDoesNotExist(self.files.botm_matching_file);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                clusteringMatched.time = toc(t_BOTMMatching);
                save(self.files.botm_matching_file, 'clusteringMatched', '-v7.3');
                disp('BOTM matching done.')
            end
            self.buffer.clusteringMatched = clusteringMatched;
        end
        
        %% --------------------------------------------------------------------
        function mergeClusters(self)
            
            try
                load(self.files.merge_clusters);
                disp('Already merged...');
            catch
                disp('Merging clusters...');
                t_merge = tic;
                clusteringMerged = [];
                % in self.buffer.clusteringMatched.ids the ids are going from 0 to nT with
                % 0 being not matched and noise, thus there is one "real" template less
                % than units (if any spike has an ID of 0 !!)
                realTemplateIdx = find(self.buffer.clusteringMatched.template_ids>0);
                nT = length(realTemplateIdx);
                if nT == 0
                    % we have only the noise cluster
                    D = [];
                    maxT = [];
                    groups = {};
                else
                    % resample only templates that do not have id==0
                    resampledTemplates = hdsort.waveforms.tResample(hdsort.waveforms.v2t(...
                        self.buffer.clusteringMatched.templates(realTemplateIdx,:), self.buffer.nC), 3, 1);
                    
                    % the indices in the groups are indices into
                    % realTemplateIdx which is index into self.buffer.clusteringMatched.template_ids
                    [groups maxT D] = hdsort.waveforms.mergeTemplates(resampledTemplates, self.noise.meanNoiseStd,...
                        'maxRelativeDistance', self.P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
                        'minCorrelation', self.P.mergeTemplates.atCorrelation);
                end
                clusteringMerged.D = D;
                clusteringMerged.maxT = maxT;
                clusteringMerged.groups = groups;
                
                spikeCuttingTf = size(self.WFS, 2) / self.buffer.nC;
                
                clusteringMerged.templates = zeros(length(groups), spikeCuttingTf*self.buffer.nC);
                clusteringMerged.ids = zeros(length(self.buffer.clusteringMatched.ids),1);
                
                % now, we need to merge the indices in ids accoring to the groups in
                % groups. But groups points into realTemplateIdx !!
                for i=1:length(groups)
                    group_idx = realTemplateIdx(clusteringMerged.groups{i});
                    group_ids = self.buffer.clusteringMatched.template_ids(group_idx);
                    idx = ismember(self.buffer.clusteringMatched.ids, group_ids);
                    % get number of spikes for each template for weighted
                    % average
                    nSpikesPerTemplate = zeros(1,length(group_ids));
                    for kk=1:length(group_ids)
                        nSpikesPerTemplate(kk) = sum(self.buffer.clusteringMatched.ids==group_ids(kk));
                    end
                    clusteringMerged.ids(idx) = i;
                    clusteringMerged.templates(i,:) = nSpikesPerTemplate*self.buffer.clusteringMatched.templates(group_idx,:)/sum(nSpikesPerTemplate);
                end
                
                self.saveIfDoesNotExist(self.files.merge_clusters);
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                clusteringMerged.time = toc(t_merge);
                save(self.files.merge_clusters, 'clusteringMerged', '-v7.3');
                fprintf('Templates after merging: %d\n', length(clusteringMerged.groups));
            end
            self.buffer.clusteringMerged = clusteringMerged;
        end
        
        
        %% --------------------------------------------------------------------
        function saveGDF(self)
            
            try
                load(self.files.gdf);
                disp('GDF already saved...');
            catch
                disp('Save gdf...');
                
                LEGindices = self.LEGindices(:);
                
                gdf = [ self.buffer.clusteringMerged.ids, ...
                    self.buffer.clusteringMatched.ts, ...
                    self.buffer.clusteringMatched.amplitudes, ...
                    LEGindices(self.buffer.clusteringMatched.amplitudeChannelIdx)];
                
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                
                save(self.files.gdf, 'gdf')
            end
            self.buffer.gdf = gdf;
        end
        
        %% --------------------------------------------------------------------
        function estimateFullTemplate(self, rawDS)
            try
                load(self.files.full_template);
                disp('Templates already estimated.')
            catch
                disp('Starting template estimation...')
                cutLeft         = self.P.templateEstimation.cutLeft;
                Tf              = self.P.templateEstimation.Tf;
                MES             = rawDS.MultiElectrode;
                elGroupIndices  = self.LEGindices;
                maxSpikes       = self.P.templateEstimation.maxSpikes;
                
                t_templateEstimation = tic;
                [wfs, nSourceSpikesPerTemplateAndChannel] = hdsort.waveforms.templateEstimation(...
                    rawDS, self.buffer.gdf, ...
                    'Tf', Tf, 'cutLeft', cutLeft, 'maxSpikes', maxSpikes);
                time = toc(t_templateEstimation);
                
                if self.STOP_ME_BECAUSE_I_AM_SLOW; return; end
                
                save(self.files.full_template, 'wfs', 'cutLeft', 'Tf', 'maxSpikes', ...
                    'MES', 'elGroupIndices', 'nSourceSpikesPerTemplateAndChannel', 'time');
                disp('Template estimation finished.')
            end
            self.buffer.fullTemplates.wfs = wfs; clear wfs;
            self.buffer.fullTemplates.cutLeft = cutLeft;
            self.buffer.fullTemplates.Tf = Tf;
            self.buffer.fullTemplates.maxSpikes = maxSpikes;
            self.buffer.fullTemplates.MultiElectrode = MES;
            self.buffer.fullTemplates.elGroupIndices = elGroupIndices;
            self.buffer.fullTemplates.nSourceSpikesPerTemplateAndChannel = nSourceSpikesPerTemplateAndChannel;
        end
        
        %% --------------------------------------------------------------------
        function saveIfDoesNotExist(self, fname)
            % This function is only needed since we have the sorter
            % parallelized and it can happen that the very same sorting runs
            % twice on the same files. In this case, the first job to finish
            % should be the one that keeps the files, they should not be
            % overwritten by the slower one.
            if exist(fname, 'file')>0
                warning(sprintf('File: %s does already exist, WILL NOT BE OVERWRITTEN - sorting ended!', fname));
                self.STOP_ME_BECAUSE_I_AM_SLOW = true;
                %self.buffer.self.STOP_ME_BECAUSE_I_AM_SLOW = self.STOP_ME_BECAUSE_I_AM_SLOW;
                % Wait a bit because of the file system issue...
                pause(10);
                return
            end
        end
        
    end
end