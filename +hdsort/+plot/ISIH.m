classdef ISIH < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        hdsort.spiketrain.
        isih
        times_ms
        threshold_line
        nUnits
        bar
        
        subhdsort.plot.
    end
    
    properties
        binSize_ms
        maxlag_ms
        refPeriod_ms
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = ISIH(hdsort.spiketrain._or_gdf, varargin)
            
            if ~iscell(hdsort.spiketrain._or_gdf) && any(size(hdsort.spiketrain._or_gdf) == 1)
                hdsort.spiketrain. = {hdsort.spiketrain._or_gdf(:)};
            elseif ~iscell(hdsort.spiketrain._or_gdf)
                hdsort.spiketrain. = mysort.hdsort.spiketrain.tGdf2cell(hdsort.spiketrain._or_gdf);
            else
                hdsort.spiketrain. = hdsort.spiketrain._or_gdf;
            end
            nUnits = numel(hdsort.spiketrain.);
            nX = ceil(sqrt(nUnits));
            nY = ceil(nUnits);
            
            P.refPeriod_ms = 1.0;
            P.binSize_ms = 0.25;
            P.maxlag_ms = 20.0;
            P.xlabel = 'ISI [ms]';
            %P.title = 'ISI refractory period violations';
            self = self@hdsort.plot.PlotInterface(P, varargin{:});
            
            if nUnits > 1
                self.subhdsort.plot. = hdsort.plot.Subhdsort.plot.([nX nY], 'showOnStartup', false);
            end
            
            self.hdsort.plot.ame = 'ISIH';
            self.hdsort.spiketrain. = hdsort.spiketrain.;
            self.nUnits = nUnits;
            
            %if self.nUnits == 1
            %    self.showAh = true;
            %end
            
            [isih times_ms] = hdsort.util.isih(self.hdsort.spiketrain., 'binSize_ms', self.binSize_ms, 'maxlag_ms', self.maxlag_ms, 'Fs', self.Fs);
            self.isih = isih;
            self.times_ms = times_ms;
            
            self.show();
        end
        
        function show_(self)
            self.setColor(self.color, 1);
            
            if self.nUnits > 1
                %show_@hdsort.plot.Subhdsort.plot.(self);
                self.subhdsort.plot..setAh(self.ah);
                self.subhdsort.plot..show();
            else
                self.useDefaultValues();
            end
            
            for ii = 1:self.nUnits
                if self.nUnits > 1
                    ah = self.subhdsort.plot..getSubhdsort.plot.andle(ii);
                else
                    ah = self.ah;
                    self.setAh(ah);
                end
                
                %b1 = bar(ah, self.times_ms, self.isih(ii, :)', 'histc');
                %set(b1, 'facecolor', self.color);
                self.bar = hdsort.plot.Bar(self.times_ms, self.isih(ii, :)', 'Style', 'histc', 'ah', ah, 'color', self.color);
                
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