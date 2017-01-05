classdef FileWrapperInterface < handle
    properties
        derivedClassName
        samplesPerSecond
        
        MultiElectrode
        fullMultiElectrode
        activeChannels
        
        memoryBufferNoiseSmad
        
        info
        name
    end
    
    methods(Abstract)
        % Need to be implemented in derived class:
        getWaveform_(self, nCut, channelIndex, varargin_)
        getData_(self, idx1, idx2)
        getNSamples_(self)
    end
    
    methods
        %------------------------------------------------------------------
        function self = FileWrapperInterface(derivedClassName, samplesPerSecond, MultiElectrode)
            
            if nargin > 2 && ~isempty(MultiElectrode)
                assert(isa(MultiElectrode, 'hdsort.filewrapper.MultiElectrode'), 'MultiElectrode must be a hdsort.filewrapper.MultiElectrode!');
                self.MultiElectrode = MultiElectrode;
            else
                self.MultiElectrode = [];
            end
            assert(isempty(samplesPerSecond) || isnumeric(samplesPerSecond), 'invalid samplerate!');
            assert(length(samplesPerSecond)<2, 'samplerate must not be a vector!');
            
            self.derivedClassName = derivedClassName;
            self.samplesPerSecond = double(samplesPerSecond);
            
            self.memoryBufferNoiseSmad = [];
            self.name = '';
        end
        
        %------------------------------------------------------------------
        function varargout = subsref(self,S)
            % This function allows you to use this object as if it was a 
            % standard matlab matrix.
            
            bIsObjectArray = length(self)~=1;
            bIsIndexAccess = length(S)==1 && strcmp(S.type, '()');
            
            if ~bIsObjectArray && bIsIndexAccess && length(S.subs)==1
                assert(S.subs{1}==1, 'indexing out of bounds!')
                varargout{1} = self;
                return
            end
            if bIsObjectArray || ~bIsIndexAccess
                % for the '.' subsref call the standard one
                [varargout{1:nargout}] = builtin('subsref', self, S);
                return
            end
            assert(strcmp(S.type, '()'), 'Only () is implemented!');
            assert(length(S.subs) > 1, '(x) is not implemented!')
            varargout{1} = self.getData(S.subs{:});
        end
        
        %------------------------------------------------------------------
        function varargout = size(self,varargin)
            dims = self.getDims();
            varargout = matlabfilecentral.parseSize.parseSize(dims,nargout,varargin{:});
        end
        
        %------------------------------------------------------------------        
        function sr = getSamplesPerSecond(self)
            sr = self.samplesPerSecond;
        end
        %------------------------------------------------------------------        
        function sr = getSampleRate(self)
            sr = self.getSamplesPerSecond();
        end
        %------------------------------------------------------------------
        function si = getSampleInterval(self)
            if isempty(self.samplesPerSecond)
                si = [];
                return
            end
            si = 1/self.samplesPerSecond;
        end
        
        %%
        %------------------------------------------------------------------
        %-------------- DataSourceInterface functions ---------------------
        %------------------------------------------------------------------
        function out = end(self, k, n)
            out = self.size(k);
        end
        %------------------------------------------------------------------
        function dims = getDims(self)
            dims = [self.getNSamples() self.getNChannels()];
        end
        
        %------------------------------------------------------------------
        function L = getNSamples(self)
            L = self.getNSamples_();
        end
        %------------------------------------------------------------------
        function n = getNChannels(self)
            n = self.MultiElectrode.getNElectrodes();
        end
        
        %------------------------------------------------------------------
        function self = transpose(self)
            error('not implemented')
        end
        %------------------------------------------------------------------
        function self = ctranspose(self)
            error('not implemented')
        end
        %------------------------------------------------------------------
        function restrictToChannels(self, channelidx)
            if isempty(self.fullMultiElectrode)
                % Make a copy of full ME before setting the sub ME as
                % active
                self.fullMultiElectrode = self.MultiElectrode;
            end
            if nargin == 1 || isempty(channelidx)
                % reset the full ME to be active
                self.MultiElectrode = self.fullMultiElectrode;
                self.activeChannels = [];
            else
                % set the active ME to be the Sub ME.
                self.MultiElectrode = self.MultiElectrode.getSubElectrode4ElIdx(channelidx);
                self.activeChannels = self.MultiElectrode.parentElectrodeIndex;
            end
        end
        %------------------------------------------------------------------
        function new = copy(self)
            tmpName = ['temp__' num2str(round(rand*10e10)) '.mat'];
            save(tmpName, 'self');
            Foo = load(tmpName);
            new = Foo.self;
            delete(tmpName);
        end

        %------------------------------------------------------------------
        function setMultiElectrode(self, ME)
            assert(ME.getNElectrodes() == self.getNChannels(), 'The multielectrode must have the same number of channels as me!');
            self.MultiElectrode = ME;
        end
        %------------------------------------------------------------------
        function ME = getMultiElectrode(self)
            ME = self.MultiElectrode;
        end
        %------------------------------------------------------------------
        function b = hasMultiElectrode(self)
            b = ~isempty(self.MultiElectrode);
        end
        
        %------------------------------------------------------------------
        function wf = getWaveform(self, t, cutLeft, cutLength, channelindex)
            assert(~isempty(t), 't must not be empty!');
            assert( isa(t, 'double'), 'Spiketimes need to be in double format');
            if nargin == 4 || isempty(channelindex)
                channelindex = 1:self.size(2);
            end
            assert(all(channelindex > 0), 'channel index out of bounds!');
            assert(all(channelindex <= self.size(2)), 'channel index out of bounds!');              
            % convert to actual channel indices !
            if ~isempty(self.activeChannels)
                channelindex = self.activeChannels(channelindex);
            end      
          
            if nargin < 6 blockwise = true; end
            
            nT = length(t);
            wf = zeros(nT, length(channelindex)*cutLength);
            
            % Sort the spiketimes and keep the unsorting index:
            [t idx_sort] = sort(t);
            [~, idx_unsort] = sort(idx_sort);
            
            tr = round(t);
            if any(tr~=t)
                warning('t contains non integer values. Subsample shifting while spike cutting is currently not supported. t will be rounded!');
            end
            t = tr;
            
            % Cut away spikes that are at the borders
            while ~isempty(t) && (t(1) <= cutLeft)
                t = t(2:end);
                idx_sort = idx_sort(2:end);
                idx_unsort = idx_unsort(2:end);
            end
            while ~isempty(t) && (t(end)-cutLeft+cutLength+1 > self.size(1))
                t(end) = [];
                idx_sort(end) = [];
                idx_unsort(end) = [];
            end
            
            if isempty(t)
                return
            end
            
            t1 = t-cutLeft;
            t2 = t-cutLeft+cutLength-1;
            nCut = length(t1);
            if nCut < nT
                warning('Could not cut all hdsort.waveforms. replacing the ones at the edges with zeros (%d from %d cut)', nCut, nT);
            end
            
            wf_t = self.getWaveform_(nCut, channelindex, cutLength, t1, t2);
            %wf_t = self.getWaveform_(t, cutLeft, cutLength, channelindex)
            
            assert(isa(wf_t, 'double'), 'Function getWaveform_ must be return a double!');
            
            % Copy cut hdsort.waveforms.back into original order, leaving those
            % that were too close to the borders zero. Might not the best
            % way to deal with this, but hey, what can we do?            
            wf(idx_sort,:) = wf_t;
        end
        
        %------------------------------------------------------------------
        % Do this in an individual function since certain
        % DataSourceInterface Implementations might want to overwrite the
        % exact way in which the file is accessed. E.g., the H5 matrix,
        % cannot access an irregular index set and has to do this either
        % chunked or individual
%         function wf = getWaveform_(self, nCut, channelindex, cutLength, t1, t2)
%             wf = zeros(nCut, length(channelindex)*cutLength);
%             % Build complete index set and access data with a single
%             % entry
%             IDX = zeros(1, nCut*cutLength);
%             for i = 1:nCut
%                 idx_idx = (i-1)*cutLength +1 : i*cutLength;
%                 IDX(idx_idx) = t1(i):t2(i); 
%             end 
%             X = self.getData_(IDX, channelindex);
%             X = reshape(X, [cutLength nCut length(channelindex)]);
%             for i = 1:nCut
%                 wf(i,:) = hdsort.waveforms.m2v(squeeze(X(:, i, :))');
%             end
%         end
        
        %------------------------------------------------------------------
        function X = getData(self, timeindex, channelindex)
            if nargin < 3
                channelindex = 1:self.size(2);
                if nargin < 2
                    timeindex = 1:self.size(1);
                end
            end
            if ischar(timeindex) && strcmp(timeindex, ':')
                timeindex = 1:self.size(1);
            end
            if ischar(channelindex) && strcmp(channelindex, ':')
                channelindex = 1:self.size(2);
            end
            if ~isempty(self.activeChannels)
                if ~islogical(channelindex)
                channelindex
                end
                channelindex = self.activeChannels(channelindex);
            end
            X = self.getData_(timeindex, channelindex);
            assert(isa(X, 'double'), 'Function getData_ must be return a double!');
        end
        
        %%
        %------------------------------------------------------------------
        %-------------- ExtendedDataSourceInterface functions -------------
        %------------------------------------------------------------------
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
        % todo
        function allspikes = getMergedSingleElectrodeDetectedSpikes(self, mergeSpikesMaxDist, varargin)
            [times pks] = self.detectSpikes(varargin{:});
            allspikes = sortrows([cell2mat(times)   cell2mat(pks)], 1);
            allspikes  = hdsort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, mergeSpikesMaxDist);
        end
        
        
    end
end
