classdef Gscatter3 < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        labels
        xData
        yData
        zData
    end
    
    properties
        zlabel
        markerArea
        markerType
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Gscatter3(xData, yData, zData, labels, varargin)
            
            P.markerArea = [];
            P.markerType = '';
            P.color = '';
            P.zlabel = '';
            
            self = self@hdsort.plot.PlotInterface(P, varargin{:})
            
            self.hdsort.plot.ame = 'Gscatter3';
            
            self.xData = xData;
            self.yData = yData;
            self.zData = zData;
            
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
             
            for g = 1:nGroups
                idx = self.labels == uLabels(g);
                self.hdsort.plot.bj(g) = scatter3(self.ah, self.xData(idx), self.yData(idx), self.zData(idx), self.markerArea, self.color(g,:), self.markerType{g});
                set(self.hdsort.plot.bj(g), 'MarkerFaceAlpha',.3, 'MarkerEdgeAlpha',.2)
            end
            
            self.useDefaultValues();
            
            zlabel(self.zlabel);
            
        end
        
    end
end