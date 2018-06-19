classdef MultiFileWrapper < hdsort.filewrapper.FileWrapperInterface
    % All files are contatenated and treated as one big recording with this object.
    
    properties
        nFiles
        fileNames
        fileWrapperList
        fileWrapperClassHandle
        
        activeFileIdx
        allFileLengths
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = MultiFileWrapper(derivedClassName, samplesPerSecond, MultiElectrode, fileNames, fileWrapperClassHandle)
            self = self@hdsort.filewrapper.FileWrapperInterface(derivedClassName, samplesPerSecond, MultiElectrode)
            
            if ~iscell(fileNames)
                fileNames = {fileNames};
            end
            self.fileNames = fileNames;
            self.nFiles = numel(self.fileNames);
            
            self.fileWrapperList = fileWrapperClassHandle(fileNames{1}, self, 1);
            self.MultiElectrode = self.fileWrapperList.MultiElectrode;
            
            for ii = 2:self.nFiles
                self.fileWrapperList(ii) = fileWrapperClassHandle(fileNames{ii}, self, ii);
                
                %% Check Multielectrodes:
                assert(  self.fileWrapperList(ii).MultiElectrode == self.MultiElectrode, 'Multielectrodes must be simmilar!')
                
            end
            self.setActiveFileIndex(1);
            
            self.name = '';
        end
        
        %------------------------------------------------------------------
        function L = getAllFileLength(self)
            if ~isempty(self.allFileLengths)
                L = self.allFileLengths;
                return
            end
            L = zeros(1, self.nFiles);
            for ii = 1:self.nFiles
                L(ii) = self.fileWrapperList(ii).size(1);
            end
            self.allFileLengths = L;
        end
        
        %------------------------------------------------------------------
        function [fileFirstIndex, fileLastIndex] = getFileIndices(self, idx)
            assert(idx>0 && idx <= self.nFiles, 'Active file index must be within range of files!')
            SL = self.getAllFileLength();
            fileOffsets = cumsum([0 SL]);
            
            fileFirstIndex = fileOffsets(idx)+1;
            fileLastIndex  = fileOffsets(idx+1);
        end
        
        %------------------------------------------------------------------
        function X = getRawData(self, timeIndex, channelIndex, fileIndex, varargin)
            if nargin < 4 || isempty(fileIndex)
                fileIndex = self.getActiveFileIndex();
            end
            if nargin < 3 || isempty(channelIndex)
                channelIndex = 1:self.size_(2);
            end
            if nargin < 2 || isempty(timeIndex)
                timeIndex = 1:self.size_(1);
            end
            X = self.fileWrapperList(fileIndex).getRawData(timeIndex, channelIndex, varargin{:});
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, timeIndex, channelIndex) 
            if nargin < 3 || strcmp(channelIndex, ':')
                channelIndex = 1:self.size(2);
            end
            
            if strcmp(timeIndex, ':')
                SL = self.getAllFileLength();
                timeIndex = 1:sum(SL);
            end
            
            X = zeros(length(timeIndex), length(channelIndex));
            
            totalLength = 0;
            for ii = 1:self.nFiles
                [fileFirstIndex, fileLastIndex] = self.getFileIndices(ii);
                
                timeIndexInThisSession = timeIndex(timeIndex >= fileFirstIndex & timeIndex <= fileLastIndex) - fileFirstIndex +1;
                
                if ~isempty(timeIndexInThisSession)
%                     assert(min(timeIndexInThisSession) > 0, 'Invalid Time Index! (<1)');
%                     assert(max(timeIndexInThisSession) <= size(self.fileWrapperList(i),1), 'Invalid Time Index! (>end)');                    
                    tmp = self.fileWrapperList(ii).getData_(timeIndexInThisSession, channelIndex, ii);
                    X(totalLength+1:totalLength+length(timeIndexInThisSession),:) = tmp;
                    totalLength = totalLength+length(timeIndexInThisSession);
                end
            end
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        %function wf = getWaveform_(self, t, cutLeft, cutLength, channels)
        function wfs = getWaveform_(self, nCut, channels, cutLength, t1, t2)
            
            if nargin < 5 || isempty(channels)
                channels = 1:self.size_(2);
            end
            
            wfs = zeros(nCut, cutLength*length(channels));

            for ii = 1:self.nFiles
                [fileFirstIndex, fileLastIndex] = self.getFileIndices(ii);
                
                inThisSessionIdx = t1>=fileFirstIndex & t2<=fileLastIndex;
                t1InThisSession = t1(inThisSessionIdx) - fileFirstIndex +1;
                t2InThisSession = t2(inThisSessionIdx) - fileFirstIndex +1;
                nInThisSession = sum(inThisSessionIdx);
                
                if nInThisSession       
                    
                    wfs_tmp = self.fileWrapperList(ii).getWaveform_(nInThisSession, ...
                        channels, cutLength, t1InThisSession, t2InThisSession);
                    
                    wfs(inThisSessionIdx, :) = wfs_tmp;
                end
            end
        end
        
        %------------------------------------------------------------------
        function [missingFrames, sessionHasMissingFrames] = getMissingFrameNumbers(self)
            nSessions = length(self.fileWrapperList);
            sessionHasMissingFrames = false(1,nSessions);
            
            missingFrames = self.fileWrapperList(1).getMissingFrameNumbers();
            sessionHasMissingFrames(1) = missingFrames(1).n > 0;
            for ii = 2:nSessions
                missingFrames(ii) = self.fileWrapperList(ii).getMissingFrameNumbers();
                sessionHasMissingFrames(ii) = missingFrames(ii).n > 0;
            end
        end
        
        %------------------------------------------------------------------
        function fileNames = getSourceFileNames(self)
            fileNames = self.fileNames;
        end
        
        %------------------------------------------------------------------
        function fileName = getFilenames(self, idx)
            fileName = self.fileNames{idx};
        end
        
        %------------------------------------------------------------------
        function setActiveFileIndex(self, idx)
            assert(idx>0 && idx <= self.nFiles, 'Active file index must be within range of files!')
            self.activeFileIdx = idx;
        end
        
        %------------------------------------------------------------------
        function idx = getActiveFileIndex(self)
            idx = self.activeFileIdx;
        end
        
        %------------------------------------------------------------------
        function N = getNSamples_(self)
            N = 0;
            for ii = 1:self.nFiles
                N = N + self.fileWrapperList(ii).getNSamples();
            end
        end
        
    end
end
