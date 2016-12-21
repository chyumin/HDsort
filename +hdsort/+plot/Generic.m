classdef Generic < myplot.PlotInterface
    properties (SetAccess=protected)
        %xData
        %yData
        %nLines
    end
    
    properties
        debug
        %yShift
        
        %plotAll
        %plotMean
        %plotMedian
        %plotStd
        %meanColor
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Generic(varargin)
            % Convert row vector to column vector:
            %if size(X, 1) == 1
            %    X = X(:);
            %end
            
            %if nargin == 1
            %    Y = X;
            %    X = (1:size(X,1))';
            %    varargin = {};
            %elseif ischar(Y)
            %    varargin = {Y, varargin{:}};
            %    Y = X;
            %    X = (1:size(X,1))';
            %end
            
            %if size(Y, 1) == 1
            %    Y = Y(:);
            %end
            %assert(size(X,1) == size(Y,1), 'X and Y must have the same length!')
            % 
            %P.yShift = 0.0;
            P.debug = false;
            %P.plotAll = true;
            %P.plotMean = false;
            %P.plotMedian = false;
            %P.plotStd = false;
            %P.meanColor = [0,0,0];
            
            self = self@myplot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Generic';
            
            %self.xData = X;
            %self.yData = Y;
            %self.nLines = size(self.yData, 2);
            
            self.show();
        end
        
        function show_(self)
            %self.setColor(self.color, self.nLines);
            
            %if self.plotAll
            %    for ii = 1:self.nLines
            %        plot(self.xData, self.yData(:,ii) + self.yShift * ii, self.LineSpec, ...
            %            'Color', self.color(ii,:), ...
            %            'LineWidth', self.LineWidth);
            %    end
            %end
            % 
            %if self.plotMean
            %    assert(self.yShift == 0.0, 'Do not plot curves with a shift and the mean!')
            %    plot(self.xData, mean(self.yData,2), '--', 'Color', self.meanColor, 'LineWidth', 2.0*self.LineWidth);
            %end
            %if self.plotMedian
            %    assert(self.yShift == 0.0, 'Do not plot curves with a shift and the median!')
            %    plot(self.xData, median(self.yData,2), '-.', 'Color', self.meanColor, 'LineWidth', 2.0*self.LineWidth);
            %end
            %if self.plotStd
            %    assert(self.yShift == 0.0, 'Do not plot curves with a shift and the std!')
            %    xP = [self.xData; flipud(self.xData)];
            %    yP = [mean(self.yData,2) + std(self.yData')'; flipud(mean(self.yData,2) - std(self.yData')')];
            %    patch(xP,yP,1,'facecolor',self.meanColor,...
            %        'edgecolor','none',...
            %        'facealpha',0.2);
            %end
            
            
            %self.useDefaultValues();
        end
        
    end
end