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
            
            %self.fileWrapperList = func2str(fileWrapperClassHandle).empty();
            self.fileWrapperList = fileWrapperClassHandle(fileNames{1}, self, 1)
            self.MultiElectrode = self.fileWrapperList.MultiElectrode;
            
            for ii = 2:self.nFiles
                self.fileWrapperList(ii) = fileWrapperClassHandle(fileNames{ii}, self, ii)
                
                %% Check Multielectrodes:
                assert(  self.fileWrapperList(ii).MultiElectrode == self.MultiElectrode, 'Multielectrodes must be simmilar!')
                
            end
            self.setActiveFileIndex(1);
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
            
            fileFirstIndex = 1;
            if idx > 1
                fileFirstIndex = fileOffsets(fileIdx-1)+1;
            end
            fileLastIndex  = fileOffsets(fileIdx);
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
        function X = getData_(self, timeIndex, channelIndex, fileIndex)
            if nargin == 4 && length(fileIndex) > 1
                X = zeros(size(timeIndex,1), size(channelIndex, 1));
                %X = self.getData_(timeIndex, channelIndex, fileIndex(1));
                for ii = 1:numel(fileIndex)
                    [fileFirstIndex, fileLastIndex] = self.getFileIndices(fileIndex(ii));
                    X(fileFirstIndex:fileLastIndex) = self.getData_(timeIndex, channelIndex, fileIndex(ii));
                    %X = [X; self.getData_(timeIndex, channelIndex, fileIndex(ii))];
                end
                return
            end
            
            if nargin < 4
                fileIndex = self.getActiveFileIndex();
                if nargin < 3
                    channelIndex = 1:self.size_(2);
                    if nargin < 2
                        timeIndex = 1:self.size_(1);
                    end
                end
            end
            
            if strcmp(timeIndex, ':')
                timeIndex = 1:self.fileWrapperList(fileIndex).size(1);
            end                    
            if strcmp(channelIndex, ':')
                channelIndex = 1:self.fileWrapperList(fileIndex).size(2);
            end
            X = self.fileWrapperList(fileIndex).getData(timeIndex, channelIndex);
        end
        
        
        %------------------------------------------------------------------
        function X = getDataConcatenated_(self, timeIndex, channelIndex) 
            error('Not tested!')
            if nargin < 3 || strcmp(channelIndex, ':')
                channelIndex = 1:self.size(2);
            end
            
            SL = self.getAllFileLength();
            if strcmp(timeIndex, ':')
                timeIndex = 1:sum(SL);
            end
            
            fileOffsets = cumsum([0 SL]);
            
            X = zeros(length(timeIndex), length(channelIndex));
            
            
            
            totalLength = 0;
            for ii = 1:self.nFiles
            %for i=1:length(fileOffsets)-1
                %sessionFirstIndex = fileOffsets(i)+1;
                %sessionLastIndex  = fileOffsets(i+1);
                
                [fileFirstIndex, fileLastIndex] = self.getFileIndices(ii);
                
                timeIndexInThisSession = timeIndex(timeIndex >= fileFirstIndex & timeIndex <= fileLastIndex) - fileFirstIndex +1;
                
                if ~isempty(timeIndexInThisSession)
%                     assert(min(timeIndexInThisSession) > 0, 'Invalid Time Index! (<1)');
%                     assert(max(timeIndexInThisSession) <= size(self.fileWrapperList(i),1), 'Invalid Time Index! (>end)');                    
                    tmp = self.getData_(timeIndexInThisSession, channelIndex, ii);
                    X(totalLength+1:totalLength+length(timeIndexInThisSession),:) = tmp;
                    totalLength = totalLength+length(timeIndexInThisSession);
                end
            end
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wf = getWaveform_(self, t, cutLeft, cutLength, channels)
            if nargin < 5 || isempty(channels)
                channelIndex = 1:self.size_(2);
            end
            
            wfs = zeros(length(t), cutLength*length(channelIndex));

            for ii = 1:self.nFiles
            %for i=1:length(fileOffsets)-1
            %    sessionFirstIndex = fileOffsets(i)+1;
            %    sessionLastIndex  = fileOffsets(i+1);
                
                [fileFirstIndex, fileLastIndex] = self.getFileIndices(ii);
                
                inThisSessionIdx = t>=fileFirstIndex & t <= fileLastIndex;
                tInThisSession = t(inThisSessionIdx) - fileFirstIndex +1;
                
                if ~isempty(tInThisSession)         
                    wfs_tmp = self.fileWrapperList(i).getWaveform(tInThisSession, cutLeft, cutLength, channelIndex);
                    wfs(inThisSessionIdx, :) = wfs_tmp;
                end
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
