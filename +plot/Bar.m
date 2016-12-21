classdef Bar < plot.PlotInterface
    properties (SetAccess=protected)
        data
        edges
        bars
        hlbars
    end
    
    properties
        highlightBars
        highlightColor
        horizontal
        Style
        GroupLabel
        barWidth
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Bar(edges, data, varargin)
            
            if nargin == 1
                data = edges;
                edges = [];
                varargin = {};
            elseif ischar(data)
                varargin = {data, varargin{:}};
                data = edges;
                edges = [];
            end
            
            P.Style = 'histc';
            P.LineSpec = '-';
            P.GroupLabel = [];
            P.barWidth = [];%0.8;
            P.highlightBars = [];
            P.horizontal = false;
            P.highlightColor = [1.0, 0.0, 0.0];
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Bar';
            self.data = data;
            self.edges = edges;
            
            if size(self.data, 1) == 1
                self.data = self.data(:);
            end
            
            if isempty(self.edges)
                self.edges = linspace(0, size(self.data,1),  size(self.data,1));
            end
            assert(size(self.edges, 2) == size(self.data, 1) , 'Number of edges must correspond to number of data points')
            
            if isempty(self.barWidth) && ~strcmp(self.Style, 'stacked')
                self.barWidth = 0.8;
            elseif isempty(self.barWidth) && strcmp(self.Style, 'stacked')
                self.barWidth = 1.0;
            end
            
            self.show();
        end
        
        function show_(self)
            
            nColors = size(self.data, 2);
            self.setColor(self.color, nColors);
            
            self.ah.NextPlot = 'replace';
            
            if ~self.horizontal
                self.bars = bar(self.ah, self.edges, self.data, self.barWidth, self.Style);
            else 
                self.bars = barh(self.ah, self.edges, self.data, self.barWidth, self.Style);
            end   
            
            for ii = 1:nColors
                set(self.bars(ii),'facecolor', self.color(ii,:), 'EdgeColor', self.color(ii,:));
                set(self.bars(ii), 'LineStyle', self.LineSpec);
            end
            
            if ~isempty(self.highlightBars)
                hold(self.ah,'on');
                
                x = self.data * 0.0;
                x(self.highlightBars) = self.data(self.highlightBars);
                
                if ~self.horizontal
                    self.hlbars = bar(self.ah, self.edges, x, self.Style);
                else
                    self.hlbars = barh(self.ah, self.edges, x, self.Style);
                end   
                
                set(self.hlbars, 'facecolor', self.highlightColor);
            end
            
            if ~isempty(self.GroupLabel)
                if size(self.GroupLabel, 1) == 1
                    self.GroupLabel = repmat(self.GroupLabel, size(self.edges,2), 1);
                end
                
                barbase = cumsum([zeros(size(self.data,1),1) self.data(:,1:end-1)],2);
                grouplabelpos = self.data/2 + barbase;
                for ii = 1:size(self.edges,2)
                    if ~self.horizontal
                        text(self.barWidth*1.2+self.edges(ii)*ones(1,size(self.data,2)), grouplabelpos(ii,:), self.GroupLabel(ii,:), 'HorizontalAlignment','center')
                    else
                        text(grouplabelpos(ii,:), self.barWidth*1.0+self.edges(ii)*ones(1,size(self.data,2)), self.GroupLabel(ii,:), 'HorizontalAlignment','center')
                    end
                end
            end
            
            if ~isempty(self.XTickLabel)
                set(self.ah, 'XTick', self.edges, 'XTickLabel', self.XTickLabel)
            else
                self.useDefaultValues();
            end
            
            
        end
        
    end
end