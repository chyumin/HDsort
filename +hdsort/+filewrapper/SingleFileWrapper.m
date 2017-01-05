classdef SingleFileWrapper < hdsort.filewrapper.FileWrapperInterface
    properties
        fileName
        
        file_idx
        parent
        
        bCovIsPrecalculated
    end
    
    methods
        %------------------------------------------------------------------
        function self = SingleFileWrapper(derivedClassName, samplesPerSecond, MultiElectrode, fileName, parent, file_idx)
            if nargin < 5
                parent = [];
                file_idx = 1;
            else
                assert(nargin == 6, 'You must provide a file_idx together with the parent object!');
            end
            
            self = self@hdsort.filewrapper.FileWrapperInterface(derivedClassName, samplesPerSecond, MultiElectrode)
            
            self.fileName = fileName;
            self.parent = parent;
            self.file_idx = file_idx;
            
            self.bCovIsPrecalculated = false;
            
            self.name = fileName;
            self.info = 'If this info message is displayed, there has been an error in the construction of this object!';
        end
        
        %------------------------------------------------------------------
        function wf = getScaledWaveform(self, nCut, channelIndex, varargin)
            wf = self.getWaveform(nCut, channelIndex, varargin{:})*self.getLSB();
        end
        
        %------------------------------------------------------------------
        function X = getScaledData(self, timeIndex, channelIndex)
            tmp = self.getData(timeIndex, channelIndex)*self.getLSB();
        end
        
        %------------------------------------------------------------------
        function lsb = getLSB(self)
            lsb = self.lsb;
        end
        
        %------------------------------------------------------------------
        function getNSamples_(self)
            error('Must be implemented in the derived class!')
        end
        function getData_(self)
            error('Must be implemented in the derived class!')
        end
        function getWaveform_(self)
            error('Must be implemented in the derived class!')
        end
        
    end
end
