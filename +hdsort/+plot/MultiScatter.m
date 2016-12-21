classdef MultiScatter < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        subhdsort.plot.
        labels
        data
    end
    
    properties
        dims
        dimLabel
        markerArea
        markerType
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = MultiScatter(data, labels, varargin)
            P.markerArea = [];
            P.markerType = '';
            
            P.dims = [];
            P.dimLabel = '';
            self = self@hdsort.plot.PlotInterface(P, varargin{:})
            self.data = data;
            self.labels = labels;
            
            if isempty(self.labels)
                self.labels = ones(size(self.data, 1), 1);
            end
            
            if isempty(self.dims)
                self.dims = 1:size(self.data,2);
            end
            assert(numel(self.dims) < 15, 'You cannot display more than 14 dimensions! Pleace specify the dims to display!')
            
            self.show();
        end
        
        function show_(self)
            %self.setColor(self.color, numel(self.labels));
            
            nD = numel(self.dims);
            self.subhdsort.plot. = hdsort.plot.Subhdsort.plot.((nD-1)*(nD-1), 'upperTriangle', 1, 'offsetY', .1, 'ah', self.ah);
            
            X = 1;
            Y = X + 1;
            for ii = 1:numel(self.subhdsort.plot..subhdsort.plot.andles)
                
                sah = self.subhdsort.plot..getSubhdsort.plot.andle(ii);
                p_ = hdsort.plot.Gscatter(self.data(:,X), self.data(:,Y), self.labels, 'ah', sah, ...
                    'markerArea', self.markerArea, 'markerType', self.markerType , 'MarkerFaceAlpha', self.MarkerFaceAlpha, 'MarkerEdgeAlpha', self.MarkerEdgeAlpha, ...
                    'xlabel', [self.dimLabel num2str(X)], 'ylabel', [self.dimLabel num2str(Y)]);
                
                self.hdsort.plot.bj = [self.hdsort.plot.bj; p_];
                
                Y = Y + 1;
                if Y > nD
                    X = X + 1;
                    Y = X + 1;
                end
            end
            
            % This prevents that the axis labels of the last subhdsort.plot.are
            % removed:
            self.setAh(self.ah);
        end
        
        % -----------------------------------------------------------------
        function setAxisAll(self, varargin)
            self.setAxis(varargin{:});
            self.subhdsort.plot..linkaxes();
            self.hdsort.plot.bj(1).setAxis(self.axis);
            self.hdsort.plot.bj(1).show();
        end
        
    end
end