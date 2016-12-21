classdef Waveforms < plot.PlotInterface
    properties (SetAccess=protected)
        waveforms
        wfLength
        nWaveforms
        subplots
    end
    
    properties
        classes
        
        nChannels
        channelWiseScaling
        gdf
        IDs
        maxDist
        restrict2Class
        plotMaxNWaveformsPerClass

        spacerY
        stacked
        plotMean
        plotMedian
        plotStd
        meanColor
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Waveforms(wavs, varargin)
            P.nChannels = 1;
            P.plotMean = false;
            P.plotMedian = false;
            P.channelWiseScaling = [];
            P.classes = []
            P.gdf = [];
            P.IDs = [];
            
            P.spacerY = 0;
            P.stacked = 1;
            
            P.maxDist = [];
            P.restrict2Class = [];
            
            P.plotMaxNWaveformsPerClass = 1000;
            self = self@plot.PlotInterface(P, varargin{:})
            
            assert(~isempty(wavs), 'Input must not be empty!')
            
            
            if ndims(wavs) == 3
                [self.wfLength self.nChannels self.nWaveforms] = size(wavs);
                wavs = mysort.wf.t2v(wavs);
            else
                self.nWaveforms = size(wavs,1);
                self.wfLength = size(wavs,2)/self.nChannels;
            end
            
            %% Todo: find better solution for this:
            if ~isempty(self.channelWiseScaling)
                self.channelWiseScaling = self.channelWiseScaling(:);
                x = repmat(self.channelWiseScaling', self.wfLength,1);
                y = repmat(x(:), 1, size(wavs,1));
                wavs = wavs .* y';
                clear x y
            end
            
            if isempty(self.classes)
                self.classes = unique(self.IDs);
            end
            
            if ~isempty(self.gdf)
                self.IDs = self.gdf(:,1);
            end
            if isempty(self.IDs)
                self.IDs = ones(self.nWaveforms,1);
            elseif length(self.IDs) == 1
                self.IDs = self.IDs*ones(self.nWaveforms,1);
            end
            
            if ~isempty(self.restrict2Class)
                idx = ismember(self.IDs, self.restrict2Class);
                wavs = wavs(idx,:);
                self.IDs = self.IDs(idx);
                if ~isempty(self.gdf)
                    self.gdf = self.gdf(idx,:);
                end
            end
            
            if ~isempty(P.gdf)
                if ~isempty(P.maxDist);
                    dt = [inf; diff(P.gdf(:,2))];
                    idx = dt>P.maxDist;
                    P.gdf(idx,:) = [];
                    wavs(idx,:) = [];
                    P.IDs(idx) = [];
                end
            end
            
            
            %%
            self.waveforms = wavs;
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            self.setColor(self.color, numel(self.classes));
            
            wavs = mysort.wf.v4plot(self.waveforms, self.nChannels);
            
            if ~self.stacked
                self.subplots = plot.Subplots([length(self.classes) 1], 'spacerY', self.spacerY, 'ah', self.ah);
                %self.subplots = mysort.plot.subplots([length(self.classes) 1], 'spacerY', self.spacerY);
                
                ahs = self.subplots.getSubplotHandle();
                
                glob_min = inf;
                glob_max = -inf;
                for ii=1:length(self.classes)
                    glob_min = min(glob_min, min(min(wavs(self.IDs==self.classes(ii),:))));
                    glob_max = max(glob_max, max(max(wavs(self.IDs==self.classes(ii),:))));
                    
                    p_ = self.plotOneClass(ahs(ii), wavs, self.classes(ii), self.color(ii, :));
                    self.plotObj = [self.plotObj; p_];
                    
                    if ii<length(self.classes)
                        ahs(ii).XTick = []; ahs(ii).XTickLabel = [];
                    else
                        ahs(ii).XTick = ((1:self.nChannels)-0.5)*(self.wfLength + 1)
                        ahs(ii).XTickLabel = util.num2str((1:self.nChannels))';
                        xlabel(ahs(ii), 'electrodes');
                    end
                end
                set(ahs, 'ylim', [glob_min glob_max]);
                linkaxes(ahs(ii), 'xy');
            else
                self.plotObj = self.plotOneClass(self.ah, wavs, self.classes, self.color)
                
                self.XTick = ((1:self.nChannels)-0.5)*(self.wfLength + 1)
                self.XTickLabel = util.num2str((1:self.nChannels))';
                self.xlabel = 'electrodes';
            end
        end
        
        % -----------------------------------------------------------------
        function [plotObj] = plotOneClass(self, ah, wavs, classes, color)
            %%
            plotObj = [];
            for ii=1:length(classes)
                idx = self.IDs == classes(ii);
                if ~isempty(self.plotMaxNWaveformsPerClass) && sum(idx) > self.plotMaxNWaveformsPerClass
                   idx_ = false & idx;
                    
                    [selectedIdx, selectedLogical] = util.samplePopulation(idx, self.plotMaxNWaveformsPerClass);
                    idx_(selectedIdx(selectedLogical)) = true;
                   idx = idx_;
                end
                
                p_ = plot(ah, wavs(idx,:)', 'color', color(ii, :), 'linewidth', self.LineWidth);
               [p_.Color] = deal([color(ii, :), self.Transparency]);
                plotObj = [plotObj; p_];
                
                set(ah, 'nextplot', 'add');
                if self.plotMean
                    plot(ah, mean(wavs(self.IDs == classes(ii),:), 1), 'color', 'k', 'linewidth', 3);
                end
                if self.plotMedian
                    plot(ah, median(wavs(self.IDs == classes(ii),:), 1), 'color', 'k', 'linewidth', 3);
                end
            end
            
            mima = [min(wavs(:)) max(wavs(:))];
            for ii=1:self.nChannels-1
                plot(ah, [ii*(self.wfLength+1) ii*(self.wfLength+1)], mima, ':', 'color', [.6 .6 .8]);
            end
            axis(ah, 'tight');
        end
        
    end
end