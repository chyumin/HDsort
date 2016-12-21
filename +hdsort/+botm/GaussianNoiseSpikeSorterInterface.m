
classdef GaussianNoiseSpikeSorterInterface  < mysort.sorters.SpikeSorterInterface
    properties (Constant = true)
        
    end    
    properties        
        Covest
    end
    
    methods (Abstract)
% %         From SpikeSorterInterface
% %         sorting = sort_(self, X)
    end    
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = GaussianNoiseSpikeSorterInterface(Covest, Tf, varargin)  
            self = self@mysort.sorters.SpikeSorterInterface();
            self.P.spikeDetector = mysort.detectors.threshold('method', 'none', 'threshFactor',4); 
            self.P.minCondNumber = 10000;
            self.P.diagonalLoading = 'DL'; % May be 'DL', 'DSL', 'none'
            self.P = hdsort.util.parseInputs(self.P, varargin);
            
            self.Tf = Tf;
%             if isa(NE, 'hdsort.util.NoiseEstimator') %, 'You must provide a NoiseEstimator as first argument!');
                self.Covest = Covest;
%             else
%                 self.NE = hdsort.util.NoiseEstimator(NE, Tf);
%             end
        end   
        
        %%% ------------------------------------------------------
        function hdsort.plot.lusterProjections(self, varargin)
%             [spikes classes] = self.getSpikeWaveforms(varargin{:});
%             mysort.hdsort.plot.clusterProjection(spikes, classes, self.templates, self.NE.getNoiseCovarianceMatrix(self.Tf));
        end

        %%% ------------------------------------------------------
        function hdsort.plot.ntraClusterPCA(self, varargin)
%             [spikes classes] = self.getSpikeWaveforms(varargin{:});
%             mysort.hdsort.plot.clusters(spikes, classes, self.NE.getNoiseCovarianceMatrix(self.Tf), 'templates', self.templates);
        end
        
        %%% ------------------------------------------------------
        function hdsort.plot.wClusterPCA(self, varargin)
%             [spikes classes] = self.getSpikeWaveforms(varargin{:});
%             
%             iU = self.NE.getPrewhiteningOperator(self.Tf);
%             pwSpikesX = spikes*iU;            
%             fetX = hdsort.util.dimReductionPCA(pwSpikesX, 4);
%             mysort.hdsort.plot.clustering(fetX, classes);
%             mysort.hdsort.plot.figureName('ClusterPCA'); 
        end        
    end
end