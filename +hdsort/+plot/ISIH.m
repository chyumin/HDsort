classdef ISIH < myplot.PlotInterface
    properties (SetAccess=protected)
        spiketrains
        isih
        times_ms
        threshold_line
        nUnits
        bar
        
        subplots
    end
    
    properties
        binSize_ms
        maxlag_ms
        refPeriod_ms
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = ISIH(spiketrains_or_gdf, varargin)
            
            if ~iscell(spiketrains_or_gdf) && any(size(spiketrains_or_gdf) == 1)
                spiketrains = {spiketrains_or_gdf(:)};
            elseif ~iscell(spiketrains_or_gdf)
                spiketrains = hdsort.spiketrain.tGdf2cell(spiketrains_or_gdf);
            else
                spiketrains = spiketrains_or_gdf;
            end
            nUnits = numel(spiketrains);
            nX = ceil(sqrt(nUnits));
            nY = ceil(nUnits);
            
            P.refPeriod_ms = 1.0;
            P.binSize_ms = 0.25;
            P.maxlag_ms = 20.0;
            P.xlabel = 'ISI [ms]';
            %P.title = 'ISI refractory period violations';
            self = self@myplot.PlotInterface(P, varargin{:});
            
            if nUnits > 1
                self.subplots = myplot.Subplots([nX nY], 'showOnStartup', false);
            end
            
            self.plotName = 'ISIH';
            self.spiketrains = spiketrains;
            self.nUnits = nUnits;
            
            %if self.nUnits == 1
            %    self.showAh = true;
            %end
            
            [isih times_ms] = util.isih(self.spiketrains, 'binSize_ms', self.binSize_ms, 'maxlag_ms', self.maxlag_ms, 'Fs', self.Fs);
            self.isih = isih;
            self.times_ms = times_ms;
            
            self.show();
        end
        
        function show_(self)
            self.setColor(self.color, 1);
            
            if self.nUnits > 1
                %show_@myplot.Subplots(self);
                self.subplots.setAh(self.ah);
                self.subplots.show();
            else
                self.useDefaultValues();
            end
            
            for ii = 1:self.nUnits
                if self.nUnits > 1
                    ah = self.subplots.getSubplotHandle(ii);
                else
                    ah = self.ah;
                    self.setAh(ah);
                end
                
                %b1 = bar(ah, self.times_ms, self.isih(ii, :)', 'histc');
                %set(b1, 'facecolor', self.color);
                self.bar = myplot.Bar(self.times_ms, self.isih(ii, :)', 'Style', 'histc', 'ah', ah, 'color', self.color);
                
                % Plot a refractory period line:
                if self.refPeriod_ms > self.binSize_ms
                    m = max(self.isih(ii, :));
                    self.threshold_line = line([self.refPeriod_ms,self.refPeriod_ms], [0.0, m], 'Color', 'r', 'LineStyle', ':','LineWidth', 2.0);
                    N = sum( self.isih(ii, self.times_ms(:) < self.refPeriod_ms) );
                    text( 1, m, sprintf('N = %d', N), 'FontSize', self.FontSize);
                end
            end
            
           
        end
        
    end
end