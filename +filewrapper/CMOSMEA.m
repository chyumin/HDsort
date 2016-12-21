classdef CMOSMEA < filewrapper.MultiSessionInterface
    properties
        P
        h5info
        fileNames
        % TODO: fileNames should not be here.
        fname
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = CMOSMEA(fnames, varargin)
            P.name = 'CMOSMEA';
            P.useFilter = 0;
            P.hpf = 500;
            P.lpf = 3000;
            P.filterOrder = 2;
            P.filterType = 'butter';
            P.filterName = '';
            P = util.parseInputs(P, varargin, 'error');
            
            if ~iscell(fnames)
                fnames = {fnames};
            end
            s_per_sec = hdf5read(fnames{1}, '/Sessions/Session0/sr');  
            assert(length(s_per_sec)==1, 'samples per second is an array!?');
            
            filterFactory_ = filewrapper.FilterWrapper(P.hpf, ...
                P.lpf, s_per_sec, P.filterOrder, P.filterType);
            
            sessionList_ = filewrapper.CMOSMEASession.empty(); 
            currentSessionCount = 0;
            nSessions = 0;
            for f=1:length(fnames)
                fname = fnames{f};
                h5info_ = hdf5info(fname);            

                % init sessions
                nSessions_ = length(h5info_.GroupHierarchy(1).Groups(1).Groups);
                nSessions = nSessions+nSessions_;
                for i=1:nSessions_
                    currentSessionCount = currentSessionCount+1;
                    sessionList_(currentSessionCount) = filewrapper.CMOSMEASession(fname, h5info_,...
                        filterFactory_, i-1, P.useFilter);
                end
            end

            self = self@filewrapper.MultiSessionInterface(P.name, s_per_sec, sessionList_);
%             self = self@filewrapper.ExtendedDataSourceInterface(P.name, s_per_sec, sessionList_(1).getMultiElectrode());
            
            for i=1:nSessions
                self.sessionList(i).parent = self;
            end
            
            %% Check Multielectrodes:
            for ii = 2:nSessions
               assert(  self.sessionList(ii).MultiElectrode == self.sessionList(1).MultiElectrode, 'Multielectrodes must be simmilar!')
            end
            
            self.fileNames = fnames;
            self.fname = P.name;
            self.h5info = h5info_;
            self.P = P;
        end
        
        %------------------------------------------------------------------
        function fileNames = getSourceFileNames(self) 
            fileNames = self.fileNames;
        end
        
        %------------------------------------------------------------------
        function X = getRawData(self, timeIndex, channelIndex, sessionIndex, varargin)
            if nargin < 4 || isempty(sessionIndex)
                sessionIndex = self.activeSessionIdx;
            end
            if nargin < 3 || isempty(channelIndex)
                channelIndex = 1:self.size_(2);
            end
            if nargin < 2 || isempty(timeIndex)
                timeIndex = 1:self.size_(1);
            end
            X = self.sessionList(sessionIndex).getRawData(timeIndex, channelIndex, varargin{:});
        end
        
        %------------------------------------------------------------------
        function x = isBinaryFile(self)
            assert(~isempty(self.sessionList), 'Session list empty');
            s = self.sessionList(1);
            x = isa( s.h5matrix_raw, 'filewrapper.binaryFileMatrix');
        end
        
        %------------------------------------------------------------------
        function W = getCutWaveforms(self, timepoints, cutleft, cutright, varargin)
            error('Depricated, use getWaveform instead');
            P.channels = [];
            P.session = self.activeSessionIdx;
            P.maxLoadSamples = 100000;
            P = util.parseInputs(P, varargin, 'error');

            if isempty(timepoints)
                W=[];
                return
            end
            
            L = cutright+cutleft+1;
            assert(L>0, 'cant cut negative epochs length!');
            nConnectedChannels = self.getNChannels(); % NConnectedChannels();
            if isempty(P.channels);
                P.channels = 1:nConnectedChannels;
            else
                assert(max(P.channels)<=nConnectedChannels, 'channel index out of bounds')
            end
            
            if size(timepoints,1) == 1
                timepoints = timepoints';
            elseif size(timepoints,2) > 1
                error('timepoints must be a vector!')
            end
            
%             % group timepoints, so that all data of one group can be loaded
%             n = length(timepoints);
%             T = sortrows([timepoints; (1:n)']);
%             assert(any(sort(timepoints)-timepoints ~=0), 'timepoints must be sorted!');
            
            t1 = max(1, min(timepoints)-cutleft);
            t2 = min(self.size_(1, P.session), max(timepoints)+cutright);
            X = self.getData(t1:t2, P.channels, P.session)';
            epochs = [timepoints-cutleft timepoints+cutright]-t1 + 1;
%             eL = mysort.epoch.length(epochs);
%             cL = [0; cumsum(eL)];
%             IDX = zeros(1, cL(end));            
%             for i=1:size(epochs,1)
%                 IDX(cL(i)+1:cL(i+1)) = epochs(i,1):epochs(i,2);
%             end
%             Xs = X(:,IDX);
            W = mysort.epoch.extractWaveform(X, epochs);
        end
  
        %------------------------------------------------------------------
        function [sFR, sessionHasMissingFrames] = getFrameNumbers(self)
            nSessions = length(self.sessionList);
            sessionHasMissingFrames = false(1,nSessions);
            FR = [];
            for i = 1:nSessions
                sFR(i) = self.sessionList(i).getFrameNumbers();
                sessionHasMissingFrames(i) = ~(length(sFR(i).missing_fns)==1);
            end
        end
        
        %------------------------------------------------------------------
        function sfnames = getSessionFilenames(self, varargin)
            sfnames = self.getSessionVar('sourceFname', varargin{:});
        end
        %------------------------------------------------------------------
        function messages = getSessionMessages(self, varargin)
            messages = self.getSessionVar('message', varargin{:});
        end 
    end
end
