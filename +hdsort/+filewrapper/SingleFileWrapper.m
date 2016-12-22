classdef SingleFileWrapper < hdsort.filewrapper.FileWrapperInterface
    properties
        fileName
        
        file_idx
        parent
        
        bCovIsPrecalculated
        
        
        %sourceFname
        %message
        %h5matrix_raw
        
        %session_str
        %size_buffer
        %CL
        %connected_channels
    end
    
    %methods(Abstract)
    % Need to be defined in derived class:
    %getWaveform_(nCut, channelIndex, varargin{:})
    %getData_(self, idx1, idx2)
    %getNSamples_(self)
    %end
    
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
            
            %self.info = ['This object is intended to wrap a single preprocessed CMOSMEA datafile. ' ...
            %             'It can be used as a standalone object, but it is recommended to use it ' ...
            %             'within hdsort.filewrapper.CMOSMEA, an object that provides the same functionality ' ...
            %             'but for multiple files at the same time.'];
            self.info = 'error!';
        end
        
        %------------------------------------------------------------------
        %function wf = getWaveform(self, nCut, channelIndex, varargin)
        %    wf = self.getWaveform_(nCut, channelIndex, varargin{:});
        %    wf = double(wf);
        %end
        
        %------------------------------------------------------------------
        function wf = getScaledWaveform(self, nCut, channelIndex, varargin)
            wf = self.getWaveform(nCut, channelIndex, varargin{:})*self.getLSB();
            %wf = double(wf)*self.getLSB();
        end
        
        %------------------------------------------------------------------
        %function X = getRawData(self, timeIndex, channelIndex)
        %    tmp = self.getData_(timeIndex, channelIndex);
        %    X = double(tmp);
        %end
        
        %------------------------------------------------------------------
        function X = getScaledData(self, timeIndex, channelIndex)
            tmp = self.getData(timeIndex, channelIndex)*self.getLSB();
            %X = double(tmp)*self.getLSB();
        end
        
        %------------------------------------------------------------------
        function lsb = getLSB(self)
            lsb = self.lsb;
        end
        
        %         %------------------------------------------------------------------
        %         function Cest = getCovest(self, varargin)
        %             if self.bCovIsPrecalculated
        %                 Cest = self.getCovestFromBufferFile(varargin{:});
        %             else
        %                 Cest = [];
        %                 % Try to get Cest from some session before us
        %                 if self.file_idx > 0 && ~isempty(self.parent)
        %                     ME1 = self.MultiElectrode;
        %                     i = 0;
        %                     while i<self.session_idx && isempty(Cest)
        %                         ME2 = self.parent.sessionList(i+1).MultiElectrode;
        %                         if length(ME1.electrodeNumbers) == length(ME2.electrodeNumbers) && ...
        %                            ~any(ME1.electrodeNumbers ~= ME2.electrodeNumbers)
        %                             Cest = self.parent.sessionList(i+1).getCovest(varargin{:});
        %                         end
        %                         i = i+1;
        %                     end
        %                 end
        %                 % If Cest could not be retrieved from someone else,
        %                 % calculate with super method
        %                 if isempty(Cest)
        %                     Cest = getCovest@hdsort.filewrapper.FileWrapperInterface(self, varargin{:});
        %                 end
        %                 CestS = Cest.toStruct();
        %                 save(self.preprocCov, 'CestS');
        %                 self.bCovIsPrecalculated = true;
        %             end
        %         end
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
