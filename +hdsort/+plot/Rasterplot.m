classdef Rasterhdsort.plot.< hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        gdf
        nUnits
    end
    
    properties
        unitIds
        Units
        timeToZero
        interval
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Rasterhdsort.plot.input, varargin)
            P.Units = [];
            P.unitIds = [];
            P.timeToZero = false;
            P.xlabel = 'time [s]';
            P.interval = [];
            self = self@hdsort.plot.PlotInterface(P, varargin{:})
            
            self.hdsort.plot.ame = 'Rasterhdsort.plot.;
            
            %% Create a gdf:
            if iscell(input)
                self.gdf = mysort.hdsort.spiketrain.toGdf(input);
            elseif isa(input,'lsa.Unit')
                self.gdf = []; self.Units = [];
                for U = input
                    self.gdf = [self.gdf; U.getGdf];
                    self.Units = [self.Units; U];
                end
            else
                % Convert row vector to column vector:
                if size(input, 1) == 1
                    input = input';
                end
                
                self.gdf = input;
                if size(self.gdf, 2) == 1
                    self.gdf(:,2) = self.gdf(:,1);
                    self.gdf(:,1) = ones(size(self.gdf,1), 1);
                end
            end
            assert( isa(self.gdf, 'double'), 'Spiketimes need to be in double format');
            if isempty(self.gdf)
                self.gdf = [0 0 0];
            end
            
            if isempty(self.unitIds)
                if isempty(self.Units)
                    self.unitIds = unique(self.gdf(:,1))';
                elseif isnumeric(self.Units)
                    self.unitIds = self.Units(:)';
                else
                    self.unitIds = [self.Units.ID];
                end
            end
            self.nUnits = length(self.unitIds);
            
            self.show(); 
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            if isempty(self.gdf)
               return;
            end
            
            if ~isempty(self.interval)
                idx = self.gdf(:,2) >= self.interval(1) & self.gdf(:,2) <= self.interval(2);
                self.gdf = self.gdf(idx,:);
                self.gdf(:,2) = self.gdf(:,2) - self.interval(1);
                
                if isempty(self.gdf)
                    return
                end
                
                self.interval = [self.gdf(1,2) self.gdf(end,2)];
                clear idx
            else
                self.interval = [min(self.gdf(:,2)) max(self.gdf(:,2))];
            end
            if self.timeToZero
                time0 = self.interval(1);
            else
                time0 = 0;
            end
            
            if isempty(self.color)
                self.color = self.vectorColor(1:self.nUnits);
            elseif self.color == -1
                self.color = zeros(self.nUnits,3);
            end
            
            if ~isempty(self.gdf)
                for i=1:self.nUnits
                    sp = self.gdf( self.gdf(:,1) == self.unitIds(i) ,2);
                    numspikes = length(sp);
                    
                    xx=ones(3*numspikes,1)*nan;
                    yy=ones(3*numspikes,1)*nan;
                    
                    gap = 0.1;
                    yy(1:3:3*numspikes)= i - 0.5*(1-gap);
                    yy(2:3:3*numspikes)= i + 0.5*(1-gap);
                    xx(1:3:3*numspikes) = (sp-time0)/self.Fs;
                    xx(2:3:3*numspikes) = (sp-time0)/self.Fs;
                    
                    hdsort.plot.xx, yy, 'Color', self.color(i,:), 'linewidth',1);
                end
            end
            
            if isempty(self.YTickLabel) && self.nUnits < 30
                self.YTick = 1:self.nUnits;
                self.YTickLabel = num2str(self.unitIds');
                self.YTickMode = 'manual'; 
                self.YTickLabelMode = 'manual';
            elseif ~isempty(self.YTickLabel)
                self.YTick = 1:size(self.YTickLabel, 1);
                self.YTickMode = 'manual';
                self.YTickLabelMode = 'manual';
            else
                self.YTick = [];%1:self.nUnits;
                self.YTickMode = 'manual';
                self.YTickLabelMode = 'manual';
            end
            if  isempty(self.xlabel)
                self.XTick = [];
                self.XTickLabel = [];
            else
                self.XTick = get(self.ah, 'XTick');
                self.XTickLabel = get(self.ah, 'XTickLabel');
            end
            
            self.ylim = [0.5 self.nUnits+0.5];
        end
    end
    
    
end