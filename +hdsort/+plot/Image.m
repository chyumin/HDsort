classdef Image < myplot.PlotInterface
    properties (SetAccess=protected)
        xData
        yData
        %nLines
        image
        ih
    end
    
    properties
        colormap
        normalize 
        %yShift
        
        %plotAll
        %plotMean
        %plotMedian
        %plotStd
        %meanColor
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Image(X, Y, IMG, varargin)
            % Convert row vector to column vector:
            %if size(X, 1) == 1
            %    X = X(:);
            %end
            
            if nargin == 1
                IMG = X;
                X = (1:size(IMG,1))';
                Y = (1:size(IMG,2))';
                varargin = {};
            elseif ischar(Y)
                varargin = {Y, varargin{:}};
                IMG = X;
                X = (1:size(IMG,1))';
                Y = (1:size(IMG,2))';
            end
            
            if size(X, 1) == 1
                X = X(:);
            end
            if size(Y, 1) == 1
                Y = Y(:);
            end
            
           % assert(size(X,1) == size(IMG,1), 'X and IMG must have the same length!')
           % assert(size(Y,1) == size(IMG,2), 'X and IMG must have the same length!')

            % 
            %P.yShift = 0.0;
            %P.plotAll = true;
            %P.plotMean = false;
            %P.plotMedian = false;
            %P.plotStd = false;
            %P.meanColor = [0,0,0];
            P.colormap = gray; 
            P.normalize = true;
            self = self@myplot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Image';
            
            self.xData = X;
            self.yData = Y;
            
            % Implement imagesc functionality:
            if ndims(IMG) == 2
                IMG = double(IMG);
                if ~self.normalize
                    self.image = IMG;
                else
                    maxI = max(IMG(:));
                    minI = min(IMG(:));
                    IMG = 1.0*(IMG - minI)/(maxI-minI);
                end
                
                %IMG = repmat(IMG, 1, 1, 3);
                
            end
            self.image = IMG;
            
            self.show();
        end
        
        function show_(self)
            %self.setColor(self.colormap);
            
            if ndims(self.image) == 2
                colormap(self.colormap);
                self.ih = image(self.xData, self.yData, self.image*64.0);
            else
                self.ih = image(self.xData, self.yData, self.image);
            end
            
            
            %set(self.ah, 'YDir', 'reverse')
            self.useDefaultValues();
        end
        
    end
end