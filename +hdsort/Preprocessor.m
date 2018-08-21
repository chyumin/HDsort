classdef Preprocessor < handle
    % This class can preprocess Raw CMOSMEA recordings, including spike
    % detection and spike waveform cutting as a preprocessing step for the
    % HD spike sorter.
    
    properties
        P
        LEGs
        DSlist
        name
        
        nC
        samplesPerSecond
        MultiElectrode
        
        folders
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Preprocessor(DSlist, LEGs, jobName, varargin)
            
            if ~iscell(DSlist)
                DSlist = {DSlist};
            end
            self.DSlist = DSlist;
            self.LEGs = LEGs;
            self.name = jobName;
            
            P.debug = false;
            P.forceFileDeletionIfExists = false;
            P.saveRawH5FileNameList = [];
            P.parfor = true;
            
            P.prefilter = 1;
            P.hpf = 300;
            P.lpf = 7000;
            P.fir_filterOrder = 110;
            P.gainmultiplier = 256;
            
            P.estimateNoise = true;
            P.deflation = 1;
            P.chunkSize = 1e6; % 1e6*1024 ==> 7.6GB
            P.minChunkSize = 20000;
            P.downSample = 1;
            P.chunkIndividualChannels = 0;
            P.subtractMeanOverAllChannels = true;
            P.restrictToTimePeriod = [];
            P.frameNoIdx = 1027:1028;
            P.sessionName = '/Sessions/Session0';
            P.h5path = []; % default is '/sig' (see below)
            P.save_as_binary = true;
            P.map = [];
            
            P.spikeDetection.method = '-';
            P.spikeDetection.thr = 5;
            P.spikeDetection.minDist = 25;
            P.spikeDetection.mergeSpikesMaxDist = 18; % do not go below 13 for HDMEA!
            
            P.spikeCutting.cutLeft = 15;
            P.spikeCutting.Tf = 65;
            
            P.noiseEstimation.maxSamples = 5 * 1e5;
            P.noiseEstimation.maxSamplesPerEpoch = 3*1e5;
            P.noiseEstimation.minDistFromSpikes = 100;
            
            P = hdsort.util.parseInputs(P, varargin);
            
            if isempty(P.h5path)
                P.h5path = '/sig';
            end
            
            self.P = P;
            
            % PRE CHECKS
            for i=1:length(self.DSlist)
                assert(isa(self.DSlist{i}, 'hdsort.file.FileWrapperInterface'), 'Must be a cell array of Data Sources!');
                if i==1
                    self.nC = size(self.DSlist{i},2);
                    self.samplesPerSecond = self.DSlist{i}.samplesPerSecond;
                    self.MultiElectrode = self.DSlist{i}.MultiElectrode;
                else
                    % TODO: This should be a check on the entire multielectrode
                    assert(self.nC == size(self.DSlist{i},2), 'All ds objects must have the same multielectrode !');
                end
            end
            
            self.LEGs = cellfun( @(x) x(:), self.LEGs, 'UniformOutput', false);
            channelIdx = unique(cat(1, self.LEGs{:}));
            assert( max(channelIdx) <= self.nC, 'The LEGs cannot contain channels larger than nC')
            assert( min(channelIdx) >= 1, 'The LEGs cannot contain channels < 1')
            
        end
        
        % -----------------------------------------------------------------
        function R = preprocess(self, OutFolder, preprocessorFile, forceExecution)
            % This function implements the main loop over all files and
            % all chunks. It loads a chunk as large as possible to still
            % fit into memory, and runs over all files concatenated in the
            % order given in DSlist.
            % It then calls, depending on what was specified when the class
            % was instantiated, the prefilter, writing of the prefiltered
            % files, noise estimation, spike detection, writing of the spike det file,
            % waveform cutting, writing of the waveform cut file
            %
            % Input:
            %   - OutFolder: location where the 'groups' file will be saved
            %   - preprocessorFile: a file that stores the important
            %   header information of the preprocessor (e.g. noiseStd of
            %   all channels); if the file already exists and can be
            %   loaded, the preprocessing will simply be skipped.
            %   - forceExecution [DEFAULT: false]: if true, run
            %   preprocessing anyway even if preprocessorFile already
            %   exists.
            if nargin < 4
                forceExecution = false;
            end
            
            try
                assert(~forceExecution, 'Run preprocessor even if file already exists.')
                R = load(preprocessorFile);
                return;
            catch
                disp('Start preprocessor...')
            end
            
            P_ = self.P;
            LEGs_ = self.LEGs;
            nC_ = self.nC;
            samplesPerSecond_ = self.samplesPerSecond;
            name_ = self.name;
            ME = self.MultiElectrode;
            
            %% Initialize parfor (if necessary):
            if P_.parfor
                parfor ii = 1:8 ii; end % make sure that the parallel pool is actually activated
                myCluster = parcluster('local');
                NumWorkers = myCluster.NumWorkers;
            else
                NumWorkers = 0;
            end
            
            %% Create Mea1kFileSaver:
            if ~isempty(P_.saveRawH5FileNameList)
                MS = {};
                deflation = 0;
                h5type = 'H5T_NATIVE_FLOAT';
                for fi=1:length(self.DSlist)
                    dims =  size(self.DSlist{fi});
                    maxDims = dims;
                    chunkDims = []; %[1000 dims(2)];
                    MS{fi} = hdsort.file.Mea1KFileSaver(OutFolder, P_.saveRawH5FileNameList{fi}, ...
                        h5type, dims, maxDims, chunkDims, deflation, P_.gainmultiplier, samplesPerSecond_, ...
                        ME.electrodePositions(:,1), ME.electrodePositions(:,2), ME.electrodeNumbers, ...
                        'save_as_binary', 1, 'forceFileDeletionIfExists', P_.forceFileDeletionIfExists);
                    
                    R.preprocessedFiles{fi} = fullfile(OutFolder, P_.saveRawH5FileNameList{fi});
                end
            end
            
            %% Initialize filter:
            filterSettings.gainmultiplier = P_.gainmultiplier;
            filterSettings.prefiltered = P_.prefilter;
            if P_.prefilter
                FIRb = hdsort.util.filter_design_fir(P_.hpf, P_.lpf, samplesPerSecond_, P_.fir_filterOrder);
                FIRb = FIRb(:);
                filterSettings.highpass = P_.hpf;
                filterSettings.lowpass = P_.lpf;
                filterSettings.downsamplefactor = 1;
                filterSettings.order = P_.fir_filterOrder;
                filterSettings.type = 'FIR fircls';
            else
                FIRb = [];
                filterSettings.highpass = 0;
                filterSettings.lowpass = 0;
                filterSettings.downsamplefactor = 1;
                filterSettings.order = 0;
                filterSettings.type = 'None';
            end
            R.FIRb = FIRb;
            
            bIsFirstChunk = true;
            nLEGs = length(LEGs_);
            nDS = length(self.DSlist);
            
            %% Prepare output folders:
            self.folders.output = OutFolder;
            if ~exist(self.folders.output, 'file')
                mkdir(self.folders.output);
            end
            groupFolders = cell(1, nLEGs);
            for legi = 1:nLEGs
                groupFolders{legi} = fullfile(self.folders.output, sprintf('group%04d', legi));
                if ~exist(groupFolders{legi}, 'file')
                    mkdir(groupFolders{legi});
                end
            end
            
            R.runtime.all = [];
            R.runtime.perDS = zeros(nDS, 1);
            R.runtime.perLEG = zeros(nLEGs, 1);
            
            nSpikesPerLEG = zeros(1, nLEGs);
            samplesPerDS = zeros(1,length(self.DSlist));
            
            %% Main loop over all files
            global_smad_per_channel = zeros(1,nC_);
            for fi = 1:nDS
                
                t_repTimesPerDS = tic;
                
                ds = self.DSlist{fi};
                
                combinedSpikeDetection = struct();
                combinedSpikeDetection.minDist = P_.spikeDetection.minDist;
                combinedSpikeDetection.threshold = P_.spikeDetection.thr;
                combinedSpikeDetection.spikesDetectedUp = cell(nC_, 1);
                combinedSpikeDetection.pks_up = cell(nC_, 1);
                combinedSpikeDetection.spikesDetectedDown = cell(nC_, 1);
                combinedSpikeDetection.pks_down = cell(nC_, 1);
                
                L = size(ds,1);
                samplesPerDS(fi) = L;
                chunker = hdsort.util.Chunker(L, ...
                    'chunkSize',    P_.chunkSize, ...
                    'chunkOverlap', P_.fir_filterOrder*2,...
                    'minChunkSize', P_.minChunkSize, ...
                    'progressDisplay', 'console');
                
                %% Loop over each chunk of each file
                while chunker.hasNextChunk()
                    
                    %repTimes = zeros(1, 7);
                    repTimes = array2table(zeros(1, 7), 'VariableNames', {...
                        'loadData', 'subtractMeanAll', 'detectSpikes', ...
                        'prepareParfor', 'saveData', 'noise', 'spikeCutting'});
                    repTimesPerLEG = zeros(nLEGs, 3);
                    
                    %% Load data:
                    t_loadData = tic;
                    [chunkOvp, chunk] = chunker.getNextChunk();
                    bIsLastChunk = ~chunker.hasNextChunk();
                    X = double(ds(chunkOvp(1):chunkOvp(2), :));
                    s1 = chunk(1)-chunkOvp(1)+1;
                    s2 = chunk(2)-chunk(1) + s1;
                    repTimes.loadData = toc(t_loadData);
                    
                    %% Subtract mean of all channels
                    t_subtractMeanAll = tic;
                    if P_.subtractMeanOverAllChannels
                        X = bsxfun(@minus, X, mean(X,2));
                    end
                    repTimes.subtractMeanAll = toc(t_subtractMeanAll);
                    
                    %% Spike detection:
                    t_detectSpikes = tic;
                    spikesDetectedDown = cell(nC_, 1);
                    pks_down = cell(nC_, 1);
                    
                    % Avoid passing entire P_ variable to parfor:
                    parfor_prefilter = P_.prefilter;
                    parfor_fir_filterOrder = P_.fir_filterOrder;
                    parfor_spikeDetection_minDist = P_.spikeDetection.minDist;
                    parfor_spikeDetection_thr = P_.spikeDetection.thr;
                    parfor_chunkSize = P_.chunkSize;
                    parfor (c = 1:nC_, NumWorkers)
                        %for c = 1:nC_
                        
                        X_ = X(:,c);
                        
                        % Subtract mean of this channel:
                        X_ = X_ - mean(X_);
                        
                        if parfor_prefilter
                            X_ = conv2(X_, FIRb, 'same');
                            % remove filter artifact at beginning
                            X_(1:parfor_fir_filterOrder, :) = 0;
                        end
                        
                        % Spike detection:
                        M = hdsort.file.DataMatrix(X_);
                        if bIsFirstChunk
                            % In the first chunk, compute the smad
                            [spikesDetectedDown(c), pks_down(c), ~, global_smad_per_channel(c)] = ...
                                M.detectSpikes('energyfun', @(x) -x, ...
                                    'minPeakDistance', parfor_spikeDetection_minDist, ...
                                    'thr', parfor_spikeDetection_thr, ...
                                    'chunkSize', parfor_chunkSize, ...
                                    'progressDisplay', 'none');
                        else
                            [spikesDetectedDown(c), pks_down(c)] = ...
                                M.detectSpikes('energyfun', @(x) -x, ...
                                'minPeakDistance', parfor_spikeDetection_minDist, ...
                                'thr', parfor_spikeDetection_thr, ...
                                'smad', global_smad_per_channel(c), ...
                                'chunkSize', parfor_chunkSize, ...
                                'progressDisplay', 'none');
                        end
                        X(:,c) = X_;
                    end
                    spikeDetection.spikesDetectedUp = cell(nC_,1);
                    spikeDetection.pks_up = cell(nC_,1);
                    spikeDetection.spikesDetectedDown = spikesDetectedDown;
                    spikeDetection.pks_down = pks_down;
                    
                    for c = 1:nC_
                        combinedSpikeDetection.spikesDetectedUp{c} = cat(1, combinedSpikeDetection.spikesDetectedUp{c}, spikeDetection.spikesDetectedUp{c});
                        combinedSpikeDetection.pks_up{c} = cat(1, combinedSpikeDetection.pks_up{c}, spikeDetection.pks_up{c});
                        combinedSpikeDetection.spikesDetectedDown{c} = cat(1, combinedSpikeDetection.spikesDetectedDown{c}, spikeDetection.spikesDetectedDown{c});
                        combinedSpikeDetection.pks_down{c} = cat(1, combinedSpikeDetection.pks_down{c}, spikeDetection.pks_down{c});
                    end
                    
                    if P_.debug
                        gdf_down = hdsort.spiketrain.toGdf(spikeDetection.spikesDetectedDown);
                        if ~isempty(gdf_down)
                            assert( all(gdf_down(:,2) > 0), '!!!')
                        end
                    end
                    repTimes.detectSpikes = toc(t_detectSpikes);
                    
                    
                    %% PREPARE PARFOR LOOP VARIABLES to avoid unnecessary broadcasting
                    t_prepareParfor = tic;
                    Tf = P_.spikeCutting.Tf;
                    cutLeft = P_.spikeCutting.cutLeft;
                    
                    LS = cell(1, nLEGs);
                    clear L
                    for legi = 1:nLEGs
                        
                        t0_leg = tic;
                        
                        L.els = LEGs_{legi};
                        L.nC = length(L.els);
                        L.allspikes = [cell2mat(spikeDetection.spikesDetectedUp(L.els))   cell2mat(spikeDetection.pks_up(L.els));
                            cell2mat(spikeDetection.spikesDetectedDown(L.els)) cell2mat(spikeDetection.pks_down(L.els))];
                        
                        L.groupFolder = groupFolders{legi};
                        L.X = X(:, L.els);
                        L.wfsFileName = fullfile(L.groupFolder, [name_ '.040spikes_cut.mat']);
                        L.noiseCovFile = fullfile(L.groupFolder, [name_ '.060cov.mat']);
                        
                        if bIsFirstChunk
                            % Check if files exist, in that case delete
                            if exist(L.wfsFileName, 'file')
                                delete(L.wfsFileName);
                            end
                            wfsFile = hdsort.file.WaveFormFileMat(L.wfsFileName, 'cutLeft', cutLeft, 'Tf', Tf, 'nC', L.nC);
                            if P_.debug
                                gdf = wfsFile.getGdf();
                                wfs = wfsFile(:,:);
                                assert(isempty(gdf), '!')
                                assert(isempty(wfs), '!')
                                clear gdf wfs
                            end
                            clear wfsFile
                        end
                        L.spikeCut = [];
                        LS{legi} = L;
                        
                        repTimesPerLEG(legi, 1) = toc(t0_leg);
                    end
                    maxDist = P_.spikeDetection.mergeSpikesMaxDist;
                    repTimes.prepareParfor = toc(t_prepareParfor);
                    
                    
                    %% Save X to Mea1kFileSaver:
                    t_saveData = tic;
                    if ~isempty(P_.saveRawH5FileNameList)
                        X = X';
                        MS{fi}.saveChunkTransposed(X(:,s1:s2));
                    end
                    clear X
                    repTimes.saveData = toc(t_saveData);
                    
                    %% Loop over each LEG for noise estimation:
                    t_noise = tic;
                    if bIsFirstChunk
                        parfor_minDistFromSpikes = P_.noiseEstimation.minDistFromSpikes;
                        parfor_maxSamplesPerEpoch = P_.noiseEstimation.maxSamplesPerEpoch;
                        parfor_maxSamples = P_.noiseEstimation.maxSamples;
                        
                        parfor (legi = 1:length(LS), NumWorkers)
                            %for legi = 1:length(LS)
                            Lpf1 = LS{legi};
                            
                            ts = hdsort.spiketrain.mergeSingleElectrodeDetectedSpikes(Lpf1.allspikes, maxDist);
                            % remove spikes that are in the overlap region
                            ts( (ts(:,1) + chunkOvp(1) < chunk(1)) | (ts(:,1) + chunkOvp(1) > chunk(2)), :) = [];
                            
                            % If this is the very first chunk
                            t1_leg = tic;
                            
                            if ~exist(Lpf1.noiseCovFile, 'file')
                                disp(['Estimating noise covariance in LEG ' num2str(legi) '...'])
                                noise = struct();
                                s1 = ts(:,1) - parfor_minDistFromSpikes;
                                s2 = ts(:,1) + parfor_minDistFromSpikes;
                                
                                noise.epochs = hdsort.epoch.flip(hdsort.epoch.merge([s1 s2]), size(Lpf1.X,1));
                                
                                Cest = hdsort.noise.Covest2(Lpf1.X, 'maxLag', Tf,...
                                    'maxSamplesPerEpoch', parfor_maxSamplesPerEpoch, ...
                                    'maxSamples', parfor_maxSamples, ...
                                    'noiseEpochs', noise.epochs, 'forceMethod', 'xcorr');
                                C_time_aligned_ = hdsort.noise.ccol2Cte(Cest.CCol, Tf);
                                noise.meanNoiseStd = sqrt(mean(diag(C_time_aligned_)));
                                noise.CestS = Cest.toStruct();
                                
                                m = matfile(Lpf1.noiseCovFile,'writable',true);
                                m.noise = noise;
                                disp('Done.')
                            else
                                disp('noiseCovFile already computed.')
                            end
                            
                            repTimesPerLEG(legi, 2) = toc(t1_leg);
                        end
                    end
                    repTimes.noise = toc(t_noise);
                    
                    %% Loop over each LEG for spike cutting:
                    t_spikeCutting = tic;
                    disp('Cutting spikes an all LEGs...')
                    parfor_debug = P_.debug;
                    parfor (legi = 1:length(LS), NumWorkers)
                        %for legi = 1:length(LS)
                        
                        t2_leg = tic;
                        if parfor_debug
                            disp(['Spike cutting: LEG ' num2str(legi) ' - NumWorkers: ' num2str(NumWorkers)])
                        end
                        
                        Lpf2 = LS{legi};
                        ts = hdsort.spiketrain.mergeSingleElectrodeDetectedSpikes(Lpf2.allspikes, maxDist);
                        
                        % Remove spikes that are in the overlap region:
                        ts( (ts(:,1) + chunkOvp(1) < chunk(1)) | (ts(:,1) + chunkOvp(1) > chunk(2)), :) = [];
                        
                        % Remove spikes that are too close to the beginning
                        % of the file:
                        ts(ts(:,1) < Tf,:) = [];
                        
                        % If this is the last chunk
                        if bIsLastChunk
                            % somehow, this can also occur when
                            % bIsLastChunk = false, so why bother...
                            % remove spikes that are too close to the end of the file
                            ts(ts(:,1) > size(Lpf2.X,1)-Tf,:) = [];
                        end
                        
                        % Do the actual spike cutting
                        nSp = size(ts,1);
                        if nSp > 0
                            cutEpochs = [ts(:,1) ts(:,1)+Tf-1]-cutLeft;
                            cutIdx = hdsort.epoch.toIdx(cutEpochs);
                            wfs = Lpf2.X(cutIdx, :);
                            wfs = reshape(wfs', [Lpf2.nC Tf nSp]);
                            wfs = permute(wfs, [2 1 3]);
                            wfs = hdsort.waveforms.t2v(wfs);
                            a = nSpikesPerLEG(legi);
                            
                            % Express spike times with respect to beginning
                            % of file:
                            ts_file = ts(:,1)+chunkOvp(1)-1;
                            
                            % Express spike times with respect to beginning of
                            % the sorting:
                            ts_sorting = ts_file + sum([0 samplesPerDS(1:fi-1)]);
                            
                            gdf = [ts_sorting*0+1, ts_sorting, ts(:,2)];
                            wfsFilepf = hdsort.file.WaveFormFileMat(Lpf2.wfsFileName, 'writable', true);
                            wfsFilepf.addWaveforms(wfs, gdf);
                            if parfor_debug
                                N = size(gdf, 1);
                                Wgdf = wfsFilepf.getGdf();
                                assert( isequal( Wgdf( (end-N+1):end, 1:2) , gdf(1:N, 1:2)), '!')
                            end
                            nSpikesPerLEG(legi) = a + nSp;
                        end
                        repTimesPerLEG(legi, 3) = toc(t2_leg);
                    end % parfor over LEGs
                    disp('Done.')
                    repTimes.spikeCutting = toc(t_spikeCutting);
                    
                    %%
                    R.runtime.all = [R.runtime.all; repTimes];
                    R.runtime.perLEG = R.runtime.perLEG + sum(repTimesPerLEG, 2);
                    
                    bIsFirstChunk = false;
                end % for loop over chunks
                
                missingFrames = ds.getMissingFrameNumbers();
                MS{fi}.saveMissingFrames(missingFrames);
                MS{fi}.saveDetectedSpikes(combinedSpikeDetection);
                MS{fi}.saveNoiseStd(global_smad_per_channel);
                MS{fi}.saveFilterSettings(filterSettings);
                
                R.runtime.perDS(fi) = toc(t_repTimesPerDS);
            end % for loop over files
            
            R.nSpikesPerLEG = nSpikesPerLEG;
            R.samplesPerDS = samplesPerDS;
            R.smadPerEl = global_smad_per_channel;
            R.groupFolders = groupFolders;
            
            R.preprocessor.P = P_;
            R.preprocessor.LEGs = LEGs_;
            R.preprocessor.nC = nC_;
            R.preprocessor.name = name_;
            R.preprocessor.samplesPerSecond = samplesPerSecond_;
            
            R.totalDuration = sum(R.runtime.perDS);
            save(preprocessorFile, '-struct', 'R');
        end % fun
        
    end % meth
end % class
