classdef Rateplot < myplot.PlotInterface
    properties (SetAccess=protected)
        rate
        nUnits
    end
    
    properties
        unitIds
        Units
        time
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Rateplot(rate, varargin)
            P.Units = [];
            P.unitIds = [];
            P.xlabel = 'time [s]';
            P.time = [];
            
            self = self@myplot.PlotInterface(P, varargin{:})
            
            self.plotName = 'Rateplot';
            
            % Convert row vector to column vector:
            if size(rate, 1) == 1
                rate = rate';
            end
            
            self.rate = rate;
            assert( isa(self.rate, 'double'), 'Rate need to be in double format');
            
            if isempty(self.unitIds)
                if isempty(self.Units)
                    self.unitIds = 1:size(self.rate, 2);
                elseif isnumeric(self.Units)
                    self.unitIds = self.Units(:)';
                else
                    self.unitIds = [self.Units.ID];
                end
            end
            self.nUnits = length(self.unitIds);
            
            if isempty(self.time)
                self.time = (1:size(self.rate, 1)) - 1;
                self.time = self.time(:) / self.Fs;
            else
                self.time = self.time(:);
            end
            assert(size(self.time, 1) == size(self.rate, 1), 'Time vector must be the same length as the rate vector!')
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            
            if isempty(self.color)
                self.color = hdsort.plot.vectorColor(1:self.nUnits);
            elseif self.color == -1
                self.color = zeros(self.nUnits,3);
            end
            
            maxInput = max(max(self.rate));
            for ii=1:self.nUnits
                plot(self.time, self.rate(:,ii) + maxInput*ii , 'Color', self.color(ii,:), 'linewidth', self.LineWidth);
            end
            
            if isempty(self.YTickLabel) & self.nUnits < 20
                self.YTick = 1:self.nUnits;
                self.YTickLabel = num2str(self.unitIds');
                self.ylabel = 'Unit IDs';                
            end
            if  isempty(self.xlabel)
                self.XTick = [];
                self.XTickLabel = [];
            else
                self.XTick = get(self.ah, 'XTick');
                self.XTickLabel = get(self.ah, 'XTickLabel');
            end
            
        end
    end
    
    
end