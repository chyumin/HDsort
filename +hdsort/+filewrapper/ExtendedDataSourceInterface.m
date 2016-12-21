classdef ExtendedDataSourceInterface < hdsort.filewrapper.DataSourceInterface
    properties
        memoryBufferNoiseSmad
        spikeTrains
    end
    
    methods(Abstract)
        getData_(self, idx1, idx2)
    end
    
    methods
        %------------------------------------------------------------------
        function self = ExtendedDataSourceInterface(varargin)
            self = self@hdsort.filewrapper.DataSourceInterface(varargin{:});
            self.memoryBufferNoiseSmad = [];
        end
        
        %------------------------------------------------------------------
        function Cest = getCovest(self, maxlag, maxsamples, maxdist, forceMethod)
            if nargin < 2 
                maxlag = 79;
            end
            if nargin < 3
                maxsamples = 150000;
            end
            if nargin < 4
                maxdist = 40;
            end
            if nargin < 5
                forceMethod = []; %'matmul', 'xcorr'
            end
            fprintf('Calculating Covest, that may take a while...\n');
            [times pks] = self.detectSpikes();
            times = double(cell2mat(times))';
            spikeEpochs = hdsort.epoch.merge([times(:)-50 times(:)+50]);
            noiseEpochs = hdsort.epoch.flip(spikeEpochs, size(self,1));
            t1 = tic;
            Cest = hdsort.noise.Covest2(self, 'maxLag', maxlag, ...
                'maxSamples', maxsamples, 'noiseEpochs', noiseEpochs,...
                'maxDist', maxdist, 'forceMethod', forceMethod);
            t2 = toc(t1);
            disp('Done.'); disp(t2);
        end
        %------------------------------------------------------------------
        function R = xcorr(self, varargin)
            P.channelIdx = 1:self.size(2);
            P.maxLen = 200000;
            P.maxLag = 100;            
            P.normalization = 'none';
            P = hdsort.util.parseInputs(P, varargin);
            
            R = xcorr(self(1:P.maxLen, P.channelIdx), P.maxLag, P.normalization);            
        end
        
        %------------------------------------------------------------------
        function [smad] = noiseStd(self, varargin)
            % Calculate channel wise noise standard deviation with the median
            % absolute deviation (MAD), invert data to ignore negative peaks
            % for that calculation
            P.channelIdx = 1:self.size(2);
            P.maxLen = 300000;
            P.thr = 4;
            P.Tf = 80;            
            P = hdsort.util.parseInputs(P, varargin);
           
            Len = self.size(1);
            fullChanIdx = P.channelIdx;
            if isempty(self.fullMultiElectrode)
                nC = self.MultiElectrode.getNElectrodes();
            else
                nC = self.fullMultiElectrode.getNElectrodes();
                if ~isempty(self.activeChannels)
                    fullChanIdx = self.activeChannels(P.channelIdx);
                end
            end
            
            if isempty(self.memoryBufferNoiseSmad)
                self.memoryBufferNoiseSmad = nan(1, nC);
            end
            notCalcIdx = isnan(self.memoryBufferNoiseSmad(fullChanIdx));
            if any(notCalcIdx)
                cidx = P.channelIdx(notCalcIdx);
                fullcidx = fullChanIdx(cidx);
                disp('Computing noise std...'); tic
                smadL = min(Len, P.maxLen);
                smad = hdsort.noise.estimateSigma(...
                        self.getData(1:smadL, cidx), P.Tf, P.thr);
                self.memoryBufferNoiseSmad(fullcidx) = smad;
                disp('Done.'); toc        
            end
            smad = self.memoryBufferNoiseSmad(fullChanIdx);
        end
        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, varargin)
            P.channelIdx = 1:self.size(2);
            P.chunkSize = 100000;
            P.thr = 3.5;
            P.energyfun = @(x) -x;
            P.minPeakDistance = ceil(self.getSamplesPerSecond/1000); % 1ms
            P.Len = [];
            P = hdsort.util.parseInputs(P, varargin);
            
            if isempty(P.Len)
                P.Len = self.size(1);
            end
            
            % get noise std
            smad = self.noiseStd('channelIdx', P.channelIdx);
            
            % Detect spikes in the beginning of the file
            disp('Detecting spikes...'); tic
            pks = cell(length(P.channelIdx),1);
            times = pks;
            for cidx = 1:length(P.channelIdx)
                c = P.channelIdx(cidx);
                pks{c,1} = [];
                times{c,1} = [];
            end
            chunker = hdsort.util.Chunker(P.Len, 'chunkSize', P.chunkSize, ...
                'progressDisplay', 'console', 'minChunkSize', 1000, 'chunkOverlap', 2*P.minPeakDistance);
            while chunker.hasNextChunk()
                [chunkOvp chunk] = chunker.getNextChunk();
                X = double(self.getData(chunkOvp(1):chunkOvp(2), P.channelIdx));
                for cidx = 1:length(P.channelIdx)
                    c = P.channelIdx(cidx);
                    [pks_, times_] = findpeaks(P.energyfun(X(:,c)), 'MINPEAKHEIGHT', smad(cidx)*P.thr,...
                                                            'MINPEAKDISTANCE', P.minPeakDistance);
                    pks_ = X(times_,c); % get the right amplitudes! (sign!)
                    pks_ = pks_(:);
                    times_ = times_(:);
                    % remove spikes that are outside this chunk
                    rmvIdx = (times_+chunkOvp(1) < chunk(1)) | (times_+chunkOvp(1) > chunk(2));
                    pks_(rmvIdx) = [];
                    times_(rmvIdx) = [];
                    
                    pks{c,1} = [pks{c}; pks_];
                    times{c,1} = [times{c}; times_+chunkOvp(1)-1];
                end
            end
            disp('Done.'); toc    
        end
        %------------------------------------------------------------------
        function allspikes = getMergedSingleElectrodeDetectedSpikes(self, mergeSpikesMaxDist, varargin)
            [times pks] = self.detectSpikes(varargin{:});
            allspikes = sortrows([cell2mat(times)   cell2mat(pks)], 1);
            allspikes  = hdsort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, mergeSpikesMaxDist);
        end        
    end
end
