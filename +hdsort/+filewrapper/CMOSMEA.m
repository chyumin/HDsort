classdef CMOSMEA < hdsort.filewrapper.MultiFileWrapper
    % This object is intended to open files that are preprocessed such that
    % the data can be used directly by the spike-sorter.
    % It takes as input a list of preprocessed recordings that were
    % recorded with the same chip-configuration (MultiElectrode) and
    % treates them as one concatenated recording.
    
    methods
        
        %------------------------------------------------------------------
        %% CONSTRUCTOR
        function self = CMOSMEA(fileNames, varargin)
            samplesPerSecond = hdf5read(fileNames{1}, '/Sessions/Session0/sr');  
            assert(length(samplesPerSecond)==1, 'Samples can not be an array!');
            
            self = self@hdsort.filewrapper.MultiFileWrapper('CMOSMEA', samplesPerSecond, [], fileNames, @hdsort.filewrapper.CMOSMEAFile);
            
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
