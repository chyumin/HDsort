classdef ErrorBars < plot.PlotInterface
    properties (SetAccess=protected)
        data
        xVals
        bars
        hlbars
    end
    
    properties
        highlightBars
        highlightColor
        Style
        barWidth
        errorbarscolor
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = ErrorBars(xVals, data, varargin)
            
            if nargin == 1
                data = xVals;
                xVals = [];
                varargin = {};
            elseif ischar(data)
                varargin = {data, varargin{:}};
                data = xVals;
                xVals = [];
            end
            
            P.Style = 'histc';
            P.LineSpec = 'o';
            
            P.barWidth = 0.8;
            P.highlightBars = [];
            P.highlightColor = [1.0, 0.0, 0.0];
            P.errorbarscolor = [0.0, 0.0, 0.0];
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Bar';
            self.data = data;
            self.xVals = xVals;
            
            if size(self.data, 1) == 1
                self.data = self.data(:);
            end
            
            if isempty(self.xVals)
                self.xVals = linspace(0, size(self.data,1),  size(self.data,1));
            end
            assert(size(self.xVals, 2) == size(self.data, 1) , 'Number of xVals must correspond to number of data points')
            
            if strcmp(self.Style, 'stacked')
                self.barWidth = 1.0;
            end
            
            self.show();
        end
        
        function show_(self)
            
            nColors = size(self.data, 2);
            self.setColor(self.color, nColors);
            
            %self.ah.NextPlot = 'replace';
            means = cellfun(@mean, self.data)
            sdev = cellfun(@std, self.data)
            
            
            plot(self.xVals, means, self.LineSpec, 'color', self.color)
%    ylim=range(c(avg-sdev, avg+sdev)),
%    pch=19, xlab="Measurements", ylab="Mean +/- SD",
%    main="Scatter plot with std.dev error bars"
%)

            errorbar(self.xVals, means, sdev, 'color', self.errorbarscolor)
            
            %self.bars = bar(self.ah, self.xVals, self.data, self.barWidth, self.Style);
            
            
            %for ii = 1:nColors
            %    set(self.bars(ii),'facecolor', self.color(ii,:), 'EdgeColor', self.color(ii,:));
            %    set(self.bars(ii), 'LineStyle', self.LineSpec); 
            %end
            
            %if ~isempty(self.highlightBars)
            %    x = self.data * 0.0; 
            %    x(self.highlightBars) = self.data(self.highlightBars);
                
            %    self.hlbars = bar(self.ah, self.xVals, x, self.Style);
            %    set(self.hlbars, 'facecolor', self.highlightColor);
            %end
            
            self.useDefaultValues();
        end
        
    end
end