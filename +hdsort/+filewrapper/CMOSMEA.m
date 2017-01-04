classdef CMOSMEAnew < hdsort.filewrapper.MultiFileWrapper
    properties
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = CMOSMEAnew(fileNames, varargin)
            
            
            %P.name = 'CMOSMEA';
            %P.useFilter = 0;
            %P.hpf = 500;
            %P.lpf = 3000;
            %P.filterOrder = 2;
            %P.filterType = 'butter';
            %P.filterName = '';
            %P = hdsort.util.parseInputs(P, varargin, 'error');
            
            %if ~iscell(fileNames)
            %    fileNames = {fileNames};
            %end
            samplesPerSecond = hdf5read(fileNames{1}, '/Sessions/Session0/sr');  
            assert(length(samplesPerSecond)==1, 'Samples can not be an array!');
            
            self = self@hdsort.filewrapper.MultiFileWrapper('CMOSMEA', samplesPerSecond, [], fileNames, @hdsort.filewrapper.CMOSMEAFile);
            %self = self@hdsort.filewrapper.FileWrapperInterface(derivedClassName, samplesPerSecond, MultiElectrode)
            
            %self.h5info = h5info_;
            %self.P = P;
            self.info = ['This object is for loading many preprocessed CMOSMEA datafiles. ' ...
                         'The files are concatenated!'];
        end
        
        %------------------------------------------------------------------
        function x = isBinaryFile(self, idx)
            if nargin < 2
                idx = self.getActiveFileIndex();
            end
            x = self.fileWrapperList(idx).isBinaryFile();
        end
        
        %------------------------------------------------------------------
        function [sFR, sessionHasMissingFrames] = getFrameNumbers(self)
            nSessions = length(self.fileWrapperList);
            sessionHasMissingFrames = false(1,nSessions);
            FR = [];
            for i = 1:nSessions
                sFR(i) = self.fileWrapperList(i).getFrameNumbers();
                sessionHasMissingFrames(i) = ~(length(sFR(i).missing_fns)==1);
            end
        end
         
    end
end
