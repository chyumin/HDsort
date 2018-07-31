classdef Unit < handle
    
    properties
        parentPopulation
        spikeTrain
        ID
        footprint
        cutLeft
        spikeAmplitudes
        meanAmplitude
        
        MultiElectrode
        nSpikes
        QC
        filePath
        myFile
        
        buffer
        
        
        
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Unit(varargin)
            
            P.footprint = [];
            P.spikeAmplitudes = [];
            P.parentPopulation = [];
            P.spikeTrain = [];
            P.ID = [];
            P.filePath = '';
            P.cutLeft = 1;
            
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            self.parentPopulation = P.parentPopulation;
            self.MultiElectrode = self.parentPopulation.MultiElectrode;
            self.ID = P.ID;
            self.spikeTrain = P.spikeTrain;
            self.footprint = P.footprint;
            self.cutLeft = P.cutLeft;
            self.spikeAmplitudes = P.spikeAmplitudes;
            self.meanAmplitude = mean(self.spikeAmplitudes);
            
            self.nSpikes = numel(self.spikeTrain);
            
            self.createQC();
            
            if ~isempty(P.filePath)
                self.filePath = P.filePath;
                self.myFile = fullfile(self.filePath, ['neuron_' num2str(self.ID) '.mat']); 
            end
        end
        
        % -----------------------------------------------------------------
        function combinedUnit = combine(self, otherUnit, newID)
            %if numel == 2 newID = self.ID; end
            N1 = self.nSpikes; N2 = otherUnit.nSpikes;
            P.footprint = (N1*self.footprint + N2*otherUnit.footprint) / (N1+N2);
            P.spikeAmplitudes = [self.spikeAmplitudes; otherUnit.spikeAmplitudes];
            P.spikeTrain = [self.spikeTrain; otherUnit.spikeTrain];
            P.ID = newID;
            P.parentPopulation = self.parentPopulation;
            P.filePath = self.filePath;
            combinedUnit = hdsort.results.Unit(P);
        end
        
        % -----------------------------------------------------------------
        function [splitUnits] = split(self, splitIdx, newIDs)
            assert( numel(splitIdx) == self.nSpikes, 'Give a destination for each spike!')
            uSplit = unique(splitIdx); N = numel(uSplit);
            assert( N == numel(newIDs), 'You must provide each split unit with a new ID!')
            
            splitUnits = hdsort.results.Unit.empty();
            for ii = 1:N
                P.footprint = self.footprint;
                P.spikeAmplitudes = self.spikeAmplitudes(splitIdx == uSplit(ii));
                P.parentPopulation = self.parentPopulation;
                P.spikeTrain = self.spikeTrain(splitIdx  == uSplit(ii));
                P.filePath = self.filePath;
                P.ID = newIDs(ii);
                splitUnits(ii) = hdsort.results.Unit(P);
            end
        end
        
        function QC = getQC(self)
            self.createQC();
            QC = self.QC;
        end
        function createQC(self, varargin)
            P.threshold = 4.8;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            self.QC.meanAmplitude = self.meanAmplitude;
            self.QC.allAmplitudes = [self.parentPopulation.Units.meanAmplitude];
            self.QC.allNSpikes = [self.parentPopulation.Units.nSpikes];
            
            
            if any(self.spikeAmplitudes)
                
            %% Estimate lost spikes:
            self.QC.AMPfit.bins = linspace(0, max(self.spikeAmplitudes), 100);
            [y_] = histc(self.spikeAmplitudes, self.QC.AMPfit.bins);
            self.QC.AMPfit.y = y_(:)';
            
            %% See whether the distribution is near the detection threshold:
            [self.QC.AMPfit.amp, maxIdx] = max(self.QC.AMPfit.y);
            mu_ = median(self.spikeAmplitudes);
            std_ = std(self.spikeAmplitudes);
            
            if std_ == 0 || P.threshold < mu_ - 2*std_
                fit = hdsort.util.gaussian1DFit(self.QC.AMPfit.bins, self.QC.AMPfit.y);
                
                self.QC.AMPfit.yfit = fit.yFit;
                self.QC.AMPfit.mu = mu_;
                self.QC.AMPfit.std = std_;
                self.QC.AMPfit.percentageLost = 0;
            else
                % The distribution is close to the detection threshold, so
                % we'll estimate the distribution differently:
                [self.QC.AMPfit.amp, maxIdx] = max(self.QC.AMPfit.y);
                self.QC.AMPfit.mu = self.QC.AMPfit.bins(maxIdx);
                self.QC.AMPfit.std = std(self.spikeAmplitudes);
                self.QC.AMPfit.yfit = self.QC.AMPfit.amp * gaussmf(self.QC.AMPfit.bins, [self.QC.AMPfit.std, self.QC.AMPfit.mu]);
            
                binIdx = self.QC.AMPfit.bins < P.threshold;
                Afit = sum(self.QC.AMPfit.yfit);
                Alost = sum(self.QC.AMPfit.yfit(binIdx));
                self.QC.AMPfit.percentageLost = 100*(Afit-Alost)/Afit;
            end
            
            self.QC.AMPfit.error = self.QC.AMPfit.yfit - self.QC.AMPfit.y;
            self.QC.AMPfit.percentageError = 100 * norm(self.QC.AMPfit.error)/(numel(self.QC.AMPfit.error)*max(self.QC.AMPfit.yfit));
            
            if 0
                %%
                figure; 
                 plot(self.QC.AMPfit.bins, self.QC.AMPfit.y); hold on;
                 plot(self.QC.AMPfit.bins, self.QC.AMPfit.yfit)
            end
                
            %% QC summary:
            self.QC.good = self.QC.AMPfit.percentageError < 5 && self.QC.AMPfit.percentageLost < 5;
            end
            
            
            self.QC.P = P;
        end
        
        % -----------------------------------------------------------------
        function loadFile(self)
            load(self.myFile)
        end
        
        % -----------------------------------------------------------------
        
        function [fh, P] = plotQCnew(self, varargin)
            P.save = 1;
            P.resave = false;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            try
                self.QC.AMPfit
            catch
                 self.createQC()
            end
            
            % starting as a Unit method --> have to find unit ID
            selfID = find(self.parentPopulation.unitIDs == self.ID);
            
            fh = figure     
            x0=10;
            y0=10;
            width=1100;
            height=800;
            set(gcf,'units','points','position',[x0,y0,width,height])
            
            samplingRate = 20000;
            
            % plot footprint
            QCsp(1) = subplot(6,3,[1 2 4 5]);
            self.plotFootprint('ah', QCsp(1));
            title('Footprints');
              
            % part about the amplitudes
            % plot spike amplitudes through the recording
            QCsp(2) = subplot(6,3,[7 8 10 11]);
            [selfRange, rasterColors, rasterNames]  = self.getRange(0);            
            rasterUnitsTrains = {self.parentPopulation.Units(selfRange).spikeTrain};
            rasterUnitsTrains = cellfun(@(x) x/samplingRate,rasterUnitsTrains,'un',0);
            rasterUnitsAmplitudes = {self.parentPopulation.Units(selfRange).spikeAmplitudes};
            
            %plot(rasterUnitsTrains{3}/20000, rasterUnitsAmplitudes{3});
            hold on
            p = cellfun(@plot, rasterUnitsTrains,rasterUnitsAmplitudes, 'UniformOutput', false);
            hold off
            title('spike amplitudes');
            hline = refline([0 4]);
            hline.Color = 'r';
            
            % plot spike amplitudes histogram
            QCsp(3) = subplot(6,3,[9 12]); hold on;
            %bins = linspace(min(self.spikeAmplitudes), max(self.spikeAmplitudes), 100);
            %[nn,xx] = histc(self.spikeAmplitudes, bins);
            %barh(bins, nn);
            barh(self.QC.AMPfit.bins, self.QC.AMPfit.y);
            
            % fit:
            plot(self.QC.AMPfit.yfit, self.QC.AMPfit.bins, 'r--')
            
            % vertical line
            hline = refline([0 4]);
            hline.Color = 'r';
            
            % [text on chart
            %barstrings = num2str(nn);
            %text(nn+0.2,bins,barstrings,'horizontalalignment','left','verticalalignment','top', ...
            %  'Fontsize',8);
           
            title('spike amplitude histogram');
            
            % plot firing rate
            QCsp(4) = subplot(6,3,[13 14]);
            a = min(self.spikeTrain);
            b = max(self.spikeTrain);
            [fr, x] = hdsort.spiketrain.toFiringRate(self.spikeTrain/20000, [], .5, [a b]/20000); 
            plot(x, fr);
            title('spike firing rate (bin .5 sec)');

            linkaxes(QCsp([2 3]), 'y')
            linkaxes(QCsp([2 4]), 'x')
            
            % plot ISIH
            QCsp(5) = subplot(6,3,[15]);
            hdsort.plot.ISIH( self.spikeTrain, 'title', 'ISIH', 'xlabel', 'ISI [ms]', 'ah', QCsp(5));

            
            % plot rasterplot
            [selfRange, rasterColors, rasterNames]  = self.getRange(5);
            rasterUnits = {self.parentPopulation.Units(selfRange).spikeTrain};
            
            QCsp(6) = subplot(6,3,[16 17 18]);
            self.plotFootprint('ah', QCsp(1));
            hdsort.plot.Rasterplot(rasterUnits, ...
                'color', rasterColors, 'title', 'Rasterplot of neighouring units', 'ah', QCsp(6));  
            set(QCsp(6), 'YTickLabel', rasterNames);
            
            if P.save
                self.filePath = hdsort.util.convertPathToOS(self.filePath);
                plotDir = fullfile(self.filePath, 'plots');
                [dir_exists,mess,messid] = mkdir( plotDir );
                figureName = fullfile(plotDir, ['u' num2str(self.ID) '_qc']);
                
                if ~exist([figureName '.fig']) P.resave = true; end
                
                if P.resave
                    hdsort.plot.PlotInterface.saveFigure(fh, figureName)
                end
            end    
        end
        
        
        function SP = plotQC(self, saveFilePath)
            if nargin == 2
                fileName = fullfile(saveFilePath, ['qcplot_unit_' num2str(self.ID)]);
                if exist([fileName '.png'])
                    return;
                end
            end
            
            selfID = find(self.parentPopulation.unitIDs == self.ID);
            SP = hdsort.plot.Subplots([2, 3], 'spacerX', 0.1, 'spacerY', 0.16, ... 
                'title', ['Unit ID: ' num2str(self.ID) ' (unit ' num2str(selfID) ')']);
            
            % ---
            SP.add(self.plotFootprint('title', 'Footprint'), 1);
            SP.add(hdsort.plot.ISIH( self.spikeTrain, 'title', 'ISIH', 'xlabel', 'ISI [ms]'), 2);
            
            QC = self.getQC(); 
            fig0 = figure; ah = axes();
            
            % ---
            hi1 = histogram(ah, self.spikeAmplitudes, 100);
            
            SP.add(hdsort.plot.Bar(hi1.BinEdges(1:end-1), hi1.Values, ...
                'title', 'Spike amplitude distribution', 'xlabel', 'spike amplitude [\sigma_{n}]'), 3);            
            
            %% ---
            
            [selfRange, rasterColors, rasterNames]  = self.getRange(7);
            rasterUnits = {self.parentPopulation.Units(selfRange).spikeTrain};
            SP.add(hdsort.plot.Rasterplot(rasterUnits, ...
                'color', rasterColors, 'title', 'Rasterplot of neighouring units'), 4);  
            set(SP.getSubplotHandle(4), 'YTickLabel', rasterNames);
            
            %% ------
            if 0
            %% ---
            hi2 = histogram(ah, QC.allAmplitudes, 100);
            idx2 = find(hi2.BinEdges <  self.QC.meanAmplitude);
            assert( ~isempty(idx2), 'Something went wrong here!')
            SP.add(hdsort.plot.Bar(hi2.BinEdges(1:end-1), hi2.Values, 'highlightBars', idx2(end), ...
                'title', 'Amplitude distribution of whole sorting', 'xlabel', 'Amplitude'), 4);
            
            %% ---
            hi3 = histogram(ah, QC.allNSpikes, 100);
            idx3 = find(hi3.BinEdges <  self.nSpikes);
            assert( ~isempty(idx3), 'Something went wrong here!')
            SP.add(hdsort.plot.Bar(hi3.BinEdges(1:end-1), hi3.Values, 'highlightBars', idx3(end), ...
                'title', 'Number of spikes distribution', 'xlabel', 'number of spikes'), 5);
            end
            
            SP.resizeFigure(1000, 700);
            
            close(fig0);
            if nargin == 2
                %fileName = fullfile(saveFilePath, ['qcplot_unit_' num2str(self.ID)]);
                SP.saveFig(fileName)
                SP.closeFigure();
            end
            
        end
        % -----------------------------------------------------------------
        function FP = plotFootprint(self, varargin)
            ME = self.MultiElectrode;
            if isempty(ME)
                warning('hdsort.results.Units has no MultiElectrode')
                ME = self.parentPopulation.MultiElectrode;
            end
            el_pos = ME.electrodePositions(ME.electrodePositions(:,1) > -1 & ME.electrodePositions(:,2) > -1, :);
            
            FP = hdsort.plot.Waveforms2D(self.footprint, el_pos, ...
                'title', num2str(self.ID), varargin{:});
        end
        
        % -----------------------------------------------------------------
        function st = getSpikeTrain(self)
            st = self.spikeTrain;
        end
        
        function st = getSpikeTrains(self)
            % Legacy function (return cell)
            st = {self.spikeTrain};
        end
        
        % -----------------------------------------------------------------
        function gdf = getGdf(self)
            if ~isfield(self.buffer, 'gdf') || isempty(self.buffer.gdf)
                self.buffer.gdf = [repmat(self.ID, self.nSpikes, 1) self.spikeTrain self.spikeAmplitudes];
            end
            gdf = self.buffer.gdf;
        end
        
        % -----------------------------------------------------------------
        function [selfRange, rasterColors, rasterNames] = getRange(self, numb)
            selfID = find(self.parentPopulation.unitIDs == self.ID);
            rRangeB = numb;
            rRangeA = numb;
            numUnits = length(self.parentPopulation.Units);
            
            if selfID < rRangeB+1
                rRangeBef = rRangeB;
                rRangeB = selfID-1;
                rRangeA = rRangeA+(rRangeBef-rRangeB);
            elseif selfID > numUnits-rRangeA
                rRangeB = rRangeB+(rRangeA-(numUnits-selfID));
                rRangeA = numUnits-selfID;
            end
            
            selfRange = (selfID-rRangeB:selfID+rRangeA);
            rmainColor = [0 0 1];
            rbgColor = [0.55 0.55 0.55];
            rasterColorBefore = repmat(rbgColor, [rRangeB 1]);
            rasterColorAfter = repmat(rbgColor, [rRangeA 1]);
            rasterColors = [rasterColorBefore; rmainColor; rasterColorAfter];
            rasterNames = [selfID-rRangeB:selfID+rRangeA];
        end
    end
    
    methods (Static)
       function self_ = loadobj(savedObject)
           self_ = savedObject;
           % Nothing to do yet
       end
    end
end