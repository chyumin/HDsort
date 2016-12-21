classdef Generic < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        %xData
        %yData
        %nLines
    end
    
    properties
        debug
        %yShift
        
        %hdsort.plot.ll
        %hdsort.plot.ean
        %hdsort.plot.edian
        %hdsort.plot.td
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
            %P.hdsort.plot.ll = true;
            %P.hdsort.plot.ean = false;
            %P.hdsort.plot.edian = false;
            %P.hdsort.plot.td = false;
            %P.meanColor = [0,0,0];
            
            self = self@hdsort.plot.PlotInterface(P, varargin{:});
            
            self.hdsort.plot.ame = 'Generic';
            
            %self.xData = X;
            %self.yData = Y;
            %self.nLines = size(self.yData, 2);
            
            self.show();
        end
        
        function show_(self)
            %self.setColor(self.color, self.nLines);
            
            %if self.hdsort.plot.ll
            %    for ii = 1:self.nLines
            %        hdsort.plot.self.xData, self.yData(:,ii) + self.yShift * ii, self.LineSpec, ...
            %            'Color', self.color(ii,:), ...
            %            'LineWidth', self.LineWidth);
            %    end
            %end
            % 
            %if self.hdsort.plot.ean
            %    assert(self.yShift == 0.0, 'Do not hdsort.plot.curves with a shift and the mean!')
            %    hdsort.plot.self.xData, mean(self.yData,2), '--', 'Color', self.meanColor, 'LineWidth', 2.0*self.LineWidth);
            %end
            %if self.hdsort.plot.edian
            %    assert(self.yShift == 0.0, 'Do not hdsort.plot.curves with a shift and the median!')
            %    hdsort.plot.self.xData, median(self.yData,2), '-.', 'Color', self.meanColor, 'LineWidth', 2.0*self.LineWidth);
            %end
            %if self.hdsort.plot.td
            %    assert(self.yShift == 0.0, 'Do not hdsort.plot.curves with a shift and the std!')
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