classdef DataMatrix < hdsort.filewrapper.FileWrapperInterface
    % Wrapper for any data matrix that can then be used by the
    % spike-sorter.
    
    properties
        M
    end
    
    methods
        
        %------------------------------------------------------------------
        %% CONSTRUCTOR
        function self = DataMatrix(M, varargin)
            P.debug = false;
            P.name = 'DataMatrix';
            P.samplesPerSecond = 20000;
            P.MultiElectrode = [];
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            self = self@hdsort.filewrapper.FileWrapperInterface('DataMatrix', P.samplesPerSecond, P.MultiElectrode)
            
            if isempty(self.MultiElectrode)
                N = size(M, 2);
                self.MultiElectrode = hdsort.filewrapper.MultiElectrode([ones(1, N); 1:N]', [1:N]);
            end
            
            self.M = M;
            self.name = P.name;
            self.info = ['Wrapper for any data matrix that can then be used by the spike-sorter.'];
        end
        
        %------------------------------------------------------------------
        function display(self)
            disp(self)
            disp(size(self.M));
            %disp(['Restriction: [' num2str(self.dims(1)) ' ' num2str(size(self.MultiElectrode.electrodePositions,1)) ']']);
            %disp(['Add. chans.: [' num2str(self.dimAddChans(1)) ' ' num2str(self.dimAddChans(2)) ']']);
        end
        
        %------------------------------------------------------------------
        % getData_(self, idx1, idx2)
        function X = getData_(self, varargin)
            X = self.M(varargin{:});
        end
        
        %------------------------------------------------------------------
        % getWaveform_(self, nCut, channelIndex, varargin_)
        function wf = getWaveform_(self, nCut, channelindex, cutLength, t1, t2)
            wf = zeros(nCut, length(channelindex)*cutLength);
            
            % Build complete index set and access data with a single
            % entry
            IDX = zeros(1, nCut*cutLength);
            for ii = 1:nCut
                idx_idx = (ii-1)*cutLength +1 : ii*cutLength;
                IDX(idx_idx) = t1(ii):t2(ii);
            end
            
            X = self.getData_(IDX, channelindex);
            X = reshape(X, [cutLength nCut length(channelindex)]);
            
            for ii = 1:nCut
                wf(ii,:) = hdsort.waveforms.m2v(squeeze(X(:, ii, :))');
            end
        end
        
        %------------------------------------------------------------------
        % getNSamples_(self)
        function L = getNSamples_(self)
            L = size(self.M, 1);
        end
        
    end
end