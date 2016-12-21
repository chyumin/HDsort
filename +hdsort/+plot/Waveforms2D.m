classdef Waveforms2D < plot.PlotInterface
    properties (SetAccess=protected)
        wfs
        med
        electrodePositions
    end
    
    properties
        IDs
        absThreshold
        maxNumberOfChannels
        plotAllTraces
        plotMeanConfidenceWithSigma
        plotElNumbers
        maxWaveforms
        channelIdx
        scaling
        wfsButtonDownFcn
        electrodeWidth
        
        plotMean
        plotMedian
        medianColor
        medianLineWidth
        
        flipud
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Waveforms2D(wfs, electrodePositions, varargin)
            P.IDs = [];
            P.absThreshold = [];
            P.maxNumberOfChannels = [];
            P.color = [.5 .5 .5];
            P.LineWidth = 1.0;
            
            P.medianColor = [.0 .0 .0];
            P.medianLineWidth = 2.0;
            
            P.plotAllTraces = 1;
            P.plotMean = false;
            P.plotMedian = false;
            P.plotMeanConfidenceWithSigma = [];
            
            P.plotElNumbers = [];
            P.maxWaveforms = 3000;
            P.channelIdx = [];
            P.scaling = [];
            P.wfsButtonDownFcn = [];
            P.electrodeWidth = 15;
            
            
            P.xlabel = '[\mum]';
            P.ylabel = '[\mum]';
            
            P.flipud = false;
            
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Waveforms2D';
            self.wfs = wfs;
            if ~self.flipud
                self.wfs = -self.wfs;
            end
            if ndims(self.wfs) == 2
                [nT nU] = size(self.wfs);
                nC = size(electrodePositions, 1);
                
                if nU == nC
                    nC = nU;
                    nU = 1;
                    self.wfs = reshape(self.wfs, [nT, nC, nU]);
                else
                    assert( mod(nT, nC) == 0, 'Number of elements must be divisible by number of channels!')
                    self.wfs = reshape(self.wfs, [nT/nC, nC, nU]);
                end
            end
            
            self.electrodePositions = electrodePositions;
            
            % prepare IDs
            if isempty(self.IDs)
                self.IDs = 1:size(self.wfs,3);
            elseif length(self.IDs) == 1
                self.IDs = ones(1, size(self.wfs,3))*self.IDs;
            end
            
            if isempty(self.maxNumberOfChannels)
                self.maxNumberOfChannels = size(self.wfs,2);
            end
            if isempty(self.scaling)
                self.scaling = 15/max(abs(self.wfs(:)));
            end
            
            if self.plotMedian || self.plotMean
                assert(self.plotMedian ~= self.plotMean, 'Only one can be plotted!!!')
                
                if self.plotMedian
                    self.med = median( self.wfs, 3);
                elseif self.plotMean
                    self.med = mean( self.wfs, 3);
                else
                    error('Should not happen!')
                end
                self.scaling = 15/max(abs(self.med(:)));
            end
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            if isempty(self.wfs)
                return
            end
            
            % wfs is tensor time x channels x items
            % electrodePositions is matrix, first column x, second column y
            assert(isempty(self.IDs) || length(self.IDs) == 1 || length(self.IDs) == size(self.wfs,3), 'Length of IDs does not match number of items in wfs')
            self.plotChannels()
            
            if self.plotMedian || self.plotMean
                self.plotMedianFcn();
            end
            
            
            if ~isempty(self.plotElNumbers)
                for ii = 1:nC
                    x = self.electrodePositions(ii,1);
                    y = self.electrodePositions(ii,2);
                    text(x + self.electrodeWidth, y, num2str(self.plotElNumbers(ii)), 'parent', self.ah);
                end
            end
            if ~isempty(self.channelIdx)
                colors = self.vectorColor(1:length(self.channelIdx));
                for ii = 1:length(self.channelIdx)
                    rectangle('Position', [ self.electrodePositions(self.channelIdx(ii), 1) ...
                        self.electrodePositions(self.channelIdx(ii), 2)-0.5*self.electrodeWidth self.electrodeWidth self.electrodeWidth ], ...
                        'EdgeColor', colors(ii, :) );
                end
            end
            
        end
        
        % -----------------------------------------------------------------
        function plotChannels(self)
            
            uIDs = unique(self.IDs);
            nU = length(uIDs);
            
            [Tf, nC, nWf] = size(self.wfs);
            
            set(self.ah, 'NextPlot', 'add');
            
            if self.plotAllTraces
                if nU > 1
                    if self.plotMedian
                        self.color = 0.75 + 0.25* mysort.plot.vectorColor(uIDs);
                    else
                        self.color = mysort.plot.vectorColor(uIDs);
                    end
                end
                
                for u = 1:nU
                    uIdx = self.IDs == uIDs(u);
                    [m mTf mchan_] =  util.max2D(abs(self.wfs(:,:,uIdx) ), 'each column');
                    
                    nC = min(self.maxNumberOfChannels, nC);
                    mchan = mchan_(1:nC);
                    
                    pIdx = nan((Tf+1)*nC,1);
                    Y = nan((Tf+1)*nC,1);
                    
                    for ii = 1:nC
                        s1 = ((Tf+1)*(ii-1)+1);
                        idx = s1:s1+Tf-1;
                        pIdx(idx) = self.electrodePositions(mchan(ii),1) + self.electrodeWidth*(0:Tf-1)/Tf;
                        Y(idx)    = self.electrodePositions(mchan(ii),2) - squeeze(self.wfs(:,mchan(ii),uIdx))*self.scaling;
                    end
                    
                    plot(self.ah, pIdx, Y, self.LineSpec, ...
                        'Color', self.color(u,:), ...
                        'LineWidth', self.LineWidth);
                    
                end
            end
        end
        % -----------------------------------------------------------------
        
        function plotMedianFcn(self)
            
            uIDs = unique(self.IDs);
            [Tf, nC, nWf] = size(self.wfs);
            
            pIdx = nan((Tf+1)*nC,1);
            Y = nan((Tf+1)*nC,nWf);
            
            for ii = 1:nC
                s1 = ((Tf+1)*(ii-1)+1);
                idx = s1:s1+Tf-1;
                pIdx(idx, 1) = self.electrodePositions(ii,1) + self.electrodeWidth*(0:Tf-1)/Tf;
                Y(idx)    = self.electrodePositions(ii,2) - squeeze(self.med(:,ii))*self.scaling;
            end
            
            plot(self.ah, pIdx, Y, self.LineSpec, ...
                'Color', self.medianColor, ...
                'LineWidth', self.medianLineWidth);
            
            if ~isempty(self.plotMeanConfidenceWithSigma)
                var_mean = self.plotMeanConfidenceWithSigma/sqrt(nWf);
                
                plot(self.ah, pIdx, Y + 3*var_mean, self.LineSpec, ...
                    'Color', self.medianColor, ...
                    'LineWidth', self.medianLineWidth);
                
                plot(self.ah, pIdx, Y - 3*var_mean, self.LineSpec, ...
                    'Color', self.medianColor, ...
                    'LineWidth', self.medianLineWidth);
                
            end
        end        
    end    
end