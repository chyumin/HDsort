classdef Gscatter < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        labels
        xData
        yData
    end
    
    properties
        hdsort.plot.roupsTogether
        defaultLegend
        markerArea
        markerType
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Gscatter(xData, yData, labels, varargin)
            
            P.markerArea = [];
            P.hdsort.plot.roupsTogether = false;
            P.defaultLegend = false;
            P.markerType = '';
            P.color = '';
            
            self = self@hdsort.plot.PlotInterface(P, varargin{:})
            
            self.hdsort.plot.ame = 'Gscatter';
            
            self.xData = xData;
            self.yData = yData;
            
            self.labels = labels;
            if isempty(self.labels)
                self.labels = ones(1, numel(xData));
            end
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            
            uLabels = unique(self.labels);
            nGroups = numel(uLabels);
            
            
            if isempty(self.color) && isempty(self.markerType)
                [colorVec markerType] = self.vectorColor(1:nGroups);
                self.color = colorVec;
                self.markerType = markerType;
                
            elseif ~isempty(self.color) && isempty(self.markerType)
                [~, markerType] = self.vectorColor(1:nGroups);
                self.color = repmat(self.color, nGroups, 1);
                self.markerType = markerType;
            elseif isempty(self.color) && ~isempty(self.markerType)
                self.color = self.vectorColor(1:nGroups);
                self.markerType = repmat({self.markerType}, 1, nGroups);
            elseif numel(self.markerType) == 1
                self.markerType = repmat({self.markerType}, 1, nGroups);
            end
            
            if self.hdsort.plot.roupsTogether
                %for g = 1:nGroups
                %idx = self.labels == uLabels(g);
                axes(self.ah)
                p_ = gscatter(self.xData, self.yData, self.labels, self.color, self.markerType{1}, self.markerArea, self.defaultLegend);
                
                for ii = 1:numel(p_)
                    [p_(ii).Color] = deal([self.color(ii, :), self.Transparency]);
                end            
                
                self.hdsort.plot.bj = [self.hdsort.plot.bj; p_];
                %end
            else
                for g = 1:nGroups
                    idx = self.labels == uLabels(g);
                    p_ = scatter(self.ah, self.xData(idx), self.yData(idx), self.markerArea, self.color(g,:), self.markerType{g});
                    self.hdsort.plot.bj = [self.hdsort.plot.bj; p_];
                end
            end
            
            self.useDefaultValues();
            
        end
        
    end
end