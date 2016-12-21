classdef Pie < plot.PlotInterface
    properties (SetAccess=protected)
        pie
    end
    
    properties
        numbers
        %legend
        %width
        explode
        labels
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Pie(numbers, varargin)
            
            P.showAh = false;
            %P.legend = [];
            P.explode = {};
            P.labels = {};
            P.width = 500;
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Pie';
            self.numbers = numbers;
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            self.setColor(self.color, numel(self.numbers) );
            
            if ~isempty(self.explode) && ~isempty(self.labels)
                self.pie = pie(self.ah, self.numbers, self.explode, self.labels);
            elseif ~isempty(self.explode)
                self.pie = pie(self.ah, self.numbers, self.explode);
            else
                self.pie = pie(self.ah, self.numbers);
            end 
            
            set(self.fh, 'colormap', self.color);
            
            self.showLegend();
            %if ~isempty(self.legend)
            %     legend(self.legend);
            %end
        end
        
        % -----------------------------------------------------------------
        function resizeFigure(self, width)
            if nargin == 2
                self.width = width;
            end
            
            if ~isvalid(self.fh)
                self.show();
            end
            set(self.fh, 'Units', 'pixels');
            set(self.fh, 'position', [100 100 self.width self.width+30]);%*11.0/12.0])
        end
    end
end