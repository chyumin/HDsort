classdef Population < handle
    % Input
    %   gdf_or_st   - EITHER: matrix with 2 columns
    %                         first column : ID of neuron
    %                         second column: spike time with respect to some global reference time point of 0
    %                     OR: cell array with on cell per neuron
    %                         containing the spike train for this
    %                         neuron. IN THIS CASE "unitIDs" MUST
    %                         BE PROVIDED !
    %
    %  footprints - (optional) foot prints for each unit. 3 matrix
    %               with (time x electrodes x units)
    %  datasource - (optional)
    
    properties
        Units
        splitmerge
        
        buffer
        name
        
        sortingInfo
        DataSource
        MultiElectrode
        
        noiseStd
        samplingRate
        info
        fileLocation
        file
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Population(varargin)
            if nargin == 1 && isstr(varargin{1})
                fileName = varargin{1};
                S = load(fileName);
                self.fromStruct(S);
                return
            elseif nargin == 1 && isstruct(varargin{1})
                self.fromStruct(varargin{1});
                return
            end
            
            P.gdf = []; % [unitID, spikeTimes, amplitude_multiple_of_noiseStd, detectionChannel]
            P.noiseStd = []; % noiseStd for each channel
            P.footprints = [];
            P.cutLeft = [];
            P.MultiElectrode = [];
            P.fileLocation = '';
            P.name = '';
            P.sortingInfo = [];
            P.samplingRate = 20000;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            assert(size(P.gdf,2) >= 4 , 'Must have at least 4 or more columns! [unitID, spiketime, amplitude, detection_channel]');
            [spikeTrains, unitIDs] = hdsort.spiketrain.fromGdf(P.gdf);
            amplitdues = P.gdf(:,3);
            detectionChannel = P.gdf(:,4);
            [nTf, nC, nUnits] = size(P.footprints);
            assert(nUnits == length(unitIDs), 'There must be one footprint per Unit!');
            assert(nC == length(P.noiseStd), 'There must be a noiseStd value for each channel!');
            assert(nC == length(P.MultiElectrode.electrodeNumbers), 'MultiElectrode not correct!');
            
            self.MultiElectrode = P.MultiElectrode;
            self.noiseStd = P.noiseStd;
            self.fileLocation = P.fileLocation;
            self.name = P.name;
            self.sortingInfo = P.sortingInfo;
            self.samplingRate = P.samplingRate;
            
            self.Units = hdsort.results.Unit.empty(0, length(unitIDs));
            for ui = 1:length(unitIDs)
                unitID_ = unitIDs(ui);
                spikeTrain_ = spikeTrains{ui};
                footprint_ = P.footprints(:,:,ui);
                amplitudes_ = amplitdues(P.gdf(:, 1) == unitID_);
                detectionChannel_ = detectionChannel(P.gdf(:, 1) == unitID_);
                
                self.Units(ui) = hdsort.results.Unit('ID', unitID_, ...
                    'spikeTrain', spikeTrain_, ...
                    'parentPopulation', self, ...
                    'footprint', footprint_, 'cutLeft', P.cutLeft, ...
                    'spikeAmplitudes', amplitudes_, ...
                    'detectionChannel', detectionChannel_, ...
                    'fileLocation', self.fileLocation);
            end
            self.getSpikeCounts();
        
        end % End of constructor
        
        
        %------------------------------------------------------------------
        function fromStruct(self, S)
            
            %% Turn S.Units into a struct array in case they are cells:
            if iscell(S.Units)
                S_Units = cat(1,S.Units{:});
                S.Units = S_Units;
            end
            
            %% First create a MultiElectrode:
            self.MultiElectrode = hdsort.file.MultiElectrode(S.MultiElectrode);
            
            %% 
            N = numel(S.Units);
            self.Units = hdsort.results.Unit.empty(0, N);
            for ii = 1:N
                S.Units(ii).parentPopulation = self;
                self.Units(ii) = hdsort.results.Unit(S.Units(ii));
            end
            
            if ~isfield(S, 'splitmerge')
                S.splitmerge = {};
            end
            
            nSM = size(S.splitmerge,1);
            for ii = 1:nSM
                U = S.splitmerge{ii,1}; clear 'U_struct'
                for jj = 1:numel(U)
                    U(jj).parentPopulation = self;
                    U_struct(jj) = hdsort.results.Unit(U(jj));
                end
                self.splitmerge{ii,1} = U_struct;
                
                U = S.splitmerge{ii,2}; clear 'U_struct'
                for jj = 1:numel(U)
                    U(jj).parentPopulation = self;
                    U_struct(jj) = hdsort.results.Unit(U(jj));
                end
                self.splitmerge{ii,2} = U_struct;
                self.splitmerge{ii,3} = S.splitmerge{ii,3};
            end
            if nSM == 0
                self.splitmerge = {};
            end
            
            if isfield(S, 'name')
                self.name = S.name;
            end
            if isfield(S, 'sortingInfo')
                self.sortingInfo = S.sortingInfo;
            end
            if isfield(S, 'noiseStd')
                self.noiseStd = S.noiseStd;
            end
            if isfield(S, 'info')
                self.info = S.info;
            end
            if isfield(S, 'fileLocation')
                try
                    self.fileLocation = hdsort.util.convertPathToOS(S.fileLocation);
                catch
                    self.fileLocation = S.fileLocation;
                end
            end
            if isfield(S, 'samplingRate')
                self.samplingRate =  S.samplingRate;
            else
                self.samplingRate = 20000;
            end
        end
        
        %------------------------------------------------------------------
        function S = toStruct(self)
            for ii = 1:numel(self.Units)
                S.Units(ii) = self.Units(ii).toStruct();
            end
            
            nSM = size(self.splitmerge,1);
            for ii = 1:nSM
                U = self.splitmerge{ii,1}; clear 'U_struct'
                for jj = 1:numel(U)
                    U_struct(jj) = U(jj).toStruct();
                end
                S.splitmerge{ii,1} = U_struct;
                
                U = self.splitmerge{ii,2}; clear 'U_struct'
                for jj = 1:numel(U)
                    U_struct(jj) = U(jj).toStruct();
                end
                S.splitmerge{ii,2} = U_struct;
                S.splitmerge{ii,3} = self.splitmerge{ii,3};
            end
            if nSM == 0
                S.splitmerge = {};
            end
            
            S.name = self.name;
            S.sortingInfo = self.sortingInfo;
            S.noiseStd = self.noiseStd;
            S.MultiElectrode = self.MultiElectrode.toStruct();
            S.fileLocation = self.fileLocation;
            S.info = self.info;
            S.samplingRate = self.samplingRate;
        end
        
        %------------------------------------------------------------------
        function save(self, fileName)
            S = self.toStruct();
            save(fileName, '-struct', 'S', '-v7.3')
        end
        
        % -----------------------------------------------------------------
        function [U, uIdx, uI] = getUnits(self, U_)
            if nargin == 1 || isempty(U_)
                U = self.Units;
            elseif isa(U_, 'hdsort.results.Unit')
                U = U_;
            elseif islogical(U_)
                U = self.Units(U_);
            elseif isnumeric(U_) && max(U_) <= numel(self.Units)
                U = self.Units(U_);
            elseif isnumeric(U_) && max(U_) > numel(self.Units)
                U = self.getUnitFromIDs(U_);
            else
                error('Not recognised input type!')
            end
            uIdx = ismember(self.Units, U);
            uI = find(uIdx);
            assert( sum(uIdx) > 0, 'Unit does not exist!')
        end
        
        % -----------------------------------------------------------------
        function gdf = getGdf(self, Units)
            if nargin == 1 %% full gdf
                if ~isfield(self.buffer, 'gdf') || isempty(self.buffer.gdf)
                    gdfC = arrayfun( @(x) x.getGdf(), self.Units, 'UniformOutput', 0);
                    self.buffer.gdf = cat(1, gdfC{:});
                end
                gdf = self.buffer.gdf;
            else
                Units = self.getUnits(Units);
                gdfC = arrayfun( @(x) x.getGdf(), Units, 'UniformOutput', 0);
                gdf = cat(1, gdfC{:});
            end        
        end
        
        % -----------------------------------------------------------------
        function st = getSpikeTrains(self, Units)
            if nargin == 1 %% full spiketrains
                if ~isfield(self.buffer, 'spikeTrains') || isempty(self.buffer.spikeTrains)
                    self.buffer.spikeTrains = {self.Units.spikeTrain};
                end
                st = self.buffer.spikeTrains;
            else
                [Units, uIdx, uI] = self.getUnits(Units);
                st = {Units.spikeTrain};
            end
        end
        
        
        % -----------------------------------------------------------------
        function [T, cutLeft] = getFootprint(self, Units)
            if nargin == 1 %% full footprint
                if ~isfield(self.buffer, 'footprint') || isempty(self.buffer.footprint)
                    self.buffer.footprint = cat(3, self.Units.footprint);
                    self.buffer.cutLeft = cat(1, self.Units.cutLeft);
                end
                T = self.buffer.footprint;
                cutLeft = self.buffer.cutLeft;
            else
                Units = self.getUnits(Units);
                T = cat(3, Units.footprint);
                cutLeft = cat(1, Units.cutLeft);
            end 
        end
        
        % -----------------------------------------------------------------
        function unitIDs = getUnitIDs(self, Units)
            if nargin == 1 %% full footprint
                if ~isfield(self.buffer, 'unitIDs') || isempty(self.buffer.unitIDs)
                    self.buffer.unitIDs = cat(1, self.Units.ID);
                end
                unitIDs = self.buffer.unitIDs;
            else
                Units = self.getUnits(Units);
                unitIDs = cat(1, Units.ID);
            end
        end
        
        function unitIDs = unitIDs(self) % legacy function
            unitIDs = self.getUnitIDs();
        end
        
        % -----------------------------------------------------------------
        function [splitUnits] = splitUnit(self, U, splitIdx)
            [U, uIdx, uI] = self.getUnits(U);
            
            %% Determine the new IDs:
            uSplit = unique(splitIdx); N = numel(uSplit);
            nLEG = floor(U.ID/1000);
            unitIDs = self.getUnitIDs();
            existingIDs = unitIDs(floor(unitIDs/1000) == nLEG);
            newIDs = (1:N) + max(existingIDs);
            assert( ~any( floor(newIDs/1000) ~= nLEG ), 'You cannnot have more than 1000 units per LEG!') 
            
            %% Create new units:
            splitUnits = split(U, splitIdx, newIDs);
            
            %% Remove old unit, add new ones and keep the change logged in 'splitmerge'
            self.splitmerge = [self.splitmerge; {U}, {splitUnits}, {'split'}];
            self.Units(uI) = [];
            self.Units = cat(2, self.Units, splitUnits);
            self.reset();
        end
        
        % -----------------------------------------------------------------
        function [mergedUnit] = mergeUnits(self, Units)
            [Units, uIdx, uI] = self.getUnits(Units);
            
            %% Determine the new ID:
            nLEG = floor(Units(1).ID/1000);
            unitIDs = self.getUnitIDs();
            existingIDs = unitIDs(floor(unitIDs/1000) == nLEG);
            newID = max(existingIDs) + 1;
            assert( floor(newID/1000) == nLEG, 'You cannnot have more than 1000 units per LEG!') 
            
            %% Create new units:
            mergedUnit = Units(1);
            for uu = 2:numel(Units)
                mergedUnit = combine(mergedUnit, Units(uu), newID);
            end
            
            %% Remove old unit, add new ones and keep the change logged in 'splitmerge'
            self.splitmerge = [self.splitmerge; {Units}, {mergedUnit}, {'merged'}];
            self.Units(uI) = [];
            self.Units = cat(2, self.Units, mergedUnit);
            self.reset();
        end
        
        function revertSplit(self, splitUnits)
            %% Check whether there are units that have been splitted or merged more than once:
            smidx = cellfun( @(x) isequal(x, splitUnits), self.splitmerge(:,2));
            assert( any(smidx), 'Error!')
            newUnits = self.splitmerge{smidx, 2};
            furtherSplitUnits = newUnits( ~ismember(newUnits, self.Units) );
            for fsU = furtherSplitUnits
                %% Do recursive revert:
                self.revertSplit(fsU);
            end
            
            %% Remove the split or merged units and reinstate the original one:
            smidx = cellfun( @(x) isequal(x, splitUnits), self.splitmerge(:,2));
            assert( any(smidx), 'Error!')
            
            newUnits = self.splitmerge{smidx, 2};
            assert( all(ismember(newUnits, self.Units)), 'At this point, all newUnits must be in self.Units!')
            
            [~, ~, existingI] = self.getUnits(newUnits);
            self.Units(existingI) = [];
            self.Units = cat(2, self.Units, self.splitmerge{smidx, 1});
            self.splitmerge(smidx, :) = [];
            self.reset();
        end
        
        function revertMerge(self, mergedUnit)
            assert( ismember(mergedUnit, self.Units), 'Merged Unit must be part of self.Units!')
            
            smidx = cellfun( @(x) isequal(x, mergedUnit), self.splitmerge(:,2));
            assert( any(smidx), 'Error!')
            oldUnits = self.splitmerge{smidx, 1};
            
            [~, ~, existingI] = self.getUnits(mergedUnit);
            self.Units(existingI) = [];
            self.Units = cat(2, self.Units, oldUnits);
            self.splitmerge(smidx, :) = [];
            self.reset();
        end
        
        function reset(self)
            [~, sidx] = sort(cat(1,self.Units.ID));
            self.Units = self.Units(sidx);
            self.buffer = [];
        end
        
        % -----------------------------------------------------------------
        function spikeCounts = getSpikeCounts(self)
            if ~isfield(self.buffer, 'spikeCounts') || isempty(self.buffer.spikeCounts)
                self.buffer.spikeCounts = arrayfun(@(x) length(x.spikeTrain), self.Units);
            end
            spikeCounts = self.buffer.spikeCounts;
        end
        
        % -----------------------------------------------------------------
        function name = getName(self)
            name = self.name;
        end
        
        % -----------------------------------------------------------------
        function U = getUnitFromIDs(self, unitIDs)
            idx = ismember(self.unitIDs, unitIDs);
            U = self.Units(idx);
        end
        
        % -----------------------------------------------------------------
        function U = getUnitsByID(self, ID)
            uIdx = ismember(self.unitIDs, ID);
            U = self.Units(uIdx);
        end
        
        % -----------------------------------------------------------------
        function U = getGoodUnits(self, conditions)
            QCs = [];
            for U_ = self.Units
                QCs = [QCs; U_.getQC];
            end
            idx = hdsort.util.thresholdFilter(conditions, QCs, true);
            U = self.Units(idx);
        end
        
        % -----------------------------------------------------------------
        function folderName = createQCPlots(self, folderName)
            if nargin == 1
                folderName = fullfile(self.fileLocation, 'plots');
                mkdir(folderName);
            end
            
            N = numel(self.Units);
            for ii = 1:N
                disp( [num2str(ii) ' of ' num2str(N) ' QC-plots created...']);
                self.Units(ii).plotQC(folderName);
            end
        end
        
        % -----------------------------------------------------------------
        function QC = getQC(self)
            if isempty(self.buffer.QC)
                self.buffer.QC.good = []
                for ii = 1:numel(self.Units)
                    U = self.Units(ii);
                    QC_ = U.getQC();
                    self.buffer.QC.good(ii) = QC_.good;
                    %QC_
                end
                
                warning('The quality control has not yet finished!')
                self.buffer.QC.nGood = sum(self.buffer.QC.good);
            end
            QC = self.buffer.QC;
        end
        
        % -----------------------------------------------------------------
        function [unit_idx, threshold] = getHighAmplitudeUnits(self, stdThreshold)
            if nargin == 1
                stdThreshold = 5;
            end
            amp = [self.Units.meanAmplitude];
            unit_idx = amp >= stdThreshold;
        end
        
        %% -----------------------------------------------------------------
        % PLOTS
        % -----------------------------------------------------------------
        function plot(self, varargin)
            
            P.Units = [];
            P.savingFolder = '';
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            [Units] = self.getUnits(P.Units);
            
            %  ############################################################
            %%  plot init
            F = hdsort.plot.Generic();
            x0=10;
            y0=10;
            width=1100;
            height=800;
            set(F.fh,'units','points','position',[x0,y0,width,height]);
            
            %% basic numbers - Units --> TEXT rather than IMG??
            % number of units
            % # of units before / after split / merge
            
            %% other basic checks
            % are bits included?
            % do bits vary? is there homogenous spread?
            % missing frames?
            
            %  ############################################################
            %%  1. Electrode position
            spah(1) = subplot(5,5,[1 2 6 7]);
            
            plot(self.MultiElectrode.electrodePositions(:,1), ...
                self.MultiElectrode.electrodePositions(:,2), 'sk', 'MarkerFaceColor','k');
            hold on;
            title_str = strcat('Configuration (', num2str(length(self.MultiElectrode.electrodePositions)), ...
                ' selected electrodes)');
            title(title_str);
            xlabel('μm')
            ylabel('μm')
            axis([0 3850 0 2100]);
            %xl = xlim;
            %yl = ylim;
            clear title;
            hold off;
            
            %  ############################################################
            %%  2. Footprint heatmap position
            % take 3 or 5 best electrode for each unit -- electrodes where footprints
            % show the largest spread
            spah(2) = subplot(5,5,[3 4 5 8 9 10]);
            set(spah(2),'box','on');
            hold on;
            title_str = sprintf('Best 3 footprints of each unit per electrode (size of bubble) \nand position of all units with mean amplitude > 6 (red)');
            title(title_str);
            xlabel('μm')
            ylabel('μm')
            clear title_str;
            
            BestFPrange = 3;
            
            bestFp = [];
            meanGoodFp = [];
            numEl = length(self.MultiElectrode.electrodePositions(:,1));
            % for each unit find 3 footprints with largest range of numbers (most
            % prominant) and add their indexes to the array
            for ii = 1:length(Units)
                fpRange = range(Units(ii).footprint);
                [val,Indeks] = maxk(fpRange, BestFPrange);
                bestFp = [bestFp Indeks];
                if Units(ii).meanAmplitude > 6
                    meanGoodFp = [meanGoodFp Indeks];
                end
            end
            clear ii;
            
            meanGoodFp = unique(meanGoodFp);
            colorMatrix = repmat([0    0.4470    0.7410],numEl,1);
            
            for num = 1:length(meanGoodFp)
                colorMatrix(meanGoodFp(num),:,:) = [0.8500    0.3250    0.0980];
            end
            clear num;
            % part for marking individual electrodes with coor not size
            % % map the array of 3 largest indexes per unit to the sum / count per
            % % index. Size = # of electrodes
            finalFPArray = accumarray(bestFp',1,[length(self.MultiElectrode.electrodePositions) 1])';
            finalFPArray = finalFPArray + 1; % otherwise you can't construct a color matrix (need indexes > 0)
            %
            % % Find the largest number of occurrences (how many units show one of their
            % % 3 "best" footprints in this electrode
            % MaxElFp = max(finalFPArray)+1;
            %
            % % Create colors for each number of occurrence. We need at least as many
            % % colors as there is the largest number of occurances
            % colors = jet(MaxElFp);
            % colors(1,:,:) = [0.85 0.85 0.85]; % making it more obvious where no best footprint is
            % colormap(colors);
            %
            % % Construct a color matrix
            % cMatrix = colors(finalFPArray, :);
            
            % Create scatter plot
            sz = 140; % size of the mark
            scatter(self.MultiElectrode.electrodePositions(:,1), ...
                self.MultiElectrode.electrodePositions(:,2), ...
                finalFPArray*6-5, colorMatrix, 'filled' );
            hold off;
            
            
            
            %  ############################################################
            %%  3. mean spike amplitudes
            
            % distribution
            % num / % of units with mean amplitudes...
            spah(3) = subplot(5,5,[11 12 13 16 17 18]);
            
            % create histogram of mean amplitudes
            histogram([Units.meanAmplitude], 100, 'BinLimits', [0,20])
            title_str = strcat('Mean amplitude of all (', num2str(length(Units)), ') units');
            title(title_str);
            xlabel('Mean signal-to-noise deviation')
            ylabel('# of units')
            clear title_str;
            % add annotation
            hold on;
            xl = xlim;
            yl = ylim;
            line([4, 4], ylim, 'LineWidth', 2, 'Color', 'red');
            text(4.2,yl(2)-10,'cutoff amplitude','Color','red');
            set(spah(3),'box','off');
            
            % additional label of # of units
            % text(xl(2)-0.25*xl(2),yl(2)-0.1*yl(2),strcat(num2str(length(RES.Units)), ' Units'),'Color','red', 'FontSize', 18);
            
            % add avg-max amplitude
            allMax = [];
            for j=1:length(Units)
                allMax = [allMax, max(Units(j).spikeAmplitudes)];
            end
            meanMax = mean(allMax);
            if meanMax < 50
                maxText = sprintf('mean (of max \namplitudes of each unit)');
                line([meanMax, meanMax], ylim, 'LineWidth', 1, 'Color', 'green', 'LineStyle', ':');
                text(meanMax + 0.02*xl(2),yl(2)-0.15*yl(2), maxText, 'Color', 'green', 'FontSize', 10);
            end
            
            % fit distribution
            % https://ch.mathworks.com/help/stats/examples/curve-fitting-and-distribution-fitting.html
            
            % % of units
            yPl = 30;
            textDiff = 4;
            splitDiff = 15;
            
            xVec = 0:2:20;
            hc = histcounts([Units.meanAmplitude], xVec);
            for ii = 1:length(hc)
                line([xVec(ii)+0.2, xVec(ii+1)-0.2], [yPl,yPl], 'LineWidth', 1, 'Color', [60/255 60/255 60/255]);
                numberPer = sprintf('%.2f',hc(ii)/length(Units)*100);
                text(xVec(ii) + 0.5,yPl+textDiff, strcat(numberPer, '%'), ...
                    'Color',[60/255 60/255 60/255], 'FontSize', 8);
                clear numberPer;
            end
            clear ii;
            
            hold off;
            
            hcAbove6 = sum(hc(4:end));
            line([xVec(4)+0.2, xVec(11)-0.2], [yPl+splitDiff,yPl+splitDiff], 'LineWidth', 1, 'Color', 'black');
            text(xVec(5)+0.5, yPl+textDiff+splitDiff, strcat(num2str(hcAbove6/length(Units)*100), '%'), ...
                'Color','black', 'FontSize', 10);
            
            %% 4. Noise Std
            %  #############################
            spah(4) = subplot(5,5,[14 15]);
            
            hold on;
            histogram([self.noiseStd], 100)
            title('Noise StD');
            xlabel('Std');
            ylabel('# of units');
            xl = xlim;
            yl = ylim;
            maxText2 = sprintf('mean noise st. dev.');
            line([mean(self.noiseStd), mean(self.noiseStd)], ylim, 'LineWidth', 1, 'Color', 'green', 'LineStyle', ':');
            text(mean(self.noiseStd) + 0.02*xl(2),yl(2)-0.15*yl(2), maxText2, 'Color', 'green', 'FontSize', 10);
            
            hold off;
            
            %  ############################################################
            %%  5. Footprint distribution
            spah(5) = subplot(5,5,[19 20]);
            
            hold on;
            
            x = [];
            for U = Units
                stdUnit = std(U.footprint);
                stdUnitSel = stdUnit > 3*mean(stdUnit);
                x = [x sum(stdUnitSel)];
            end
            
            histogram(x, 100)
            title_str = sprintf('Footprint size (# of footprints that are \nlarger than 2*std(footprint of the unit)');
            title(title_str);
            xlabel('footprints');
            ylabel('# of units');
            %xl = xlim;
            %yl = ylim;
            % maxText2 = sprintf('mean noise st. dev.');
            % line([mean(self.noiseStd), mean(self.noiseStd)], ylim, 'LineWidth', 1, 'Color', 'green', 'LineStyle', ':');
            % text(mean(self.noiseStd) + 0.02*xl(2),yl(2)-0.15*yl(2), maxText2, 'Color', 'green', 'FontSize', 10);
            
            clear x; clear title_str;
            hold off;
            
            %  ############################################################
            %%  6. Comulative spike plot
            spah(6) = subplot(5,5,[21 22 23 24]);
            
            ST = self.getSpikeTrains(Units);
            nBins = 10000;
            srate = 20000;
            
            window = 0.1; % seconds
            bin_width = window*self.samplingRate;
            [binnedST, bin_edges, bin_centers] = hdsort.spiketrain.binSpikeTrains(ST(:), window, self.samplingRate);
            bin_centers = bin_centers - bin_centers(1);
            freq = sum(binnedST)/(bin_width/srate);
            
            bar(bin_centers, freq, 'barWidth',1);
            title('Total number of spikes (for all units)');
            xlabel('time [s]')
            ylabel('# of spikes')
            
            clear bin_edges; clear bin_width; clear bin_centers;
            clear n; clear nBins; clear srate; clear freq;
            
            %  ############################################################
            %%  7. Distribution of spike num per units
            % get total duration of the recordings:
            spah(7) = subplot(5,5,[25]);
            totalDur = 0; % in sec
            for td=1:length(self.sortingInfo.fileInfo)
                totalDur = totalDur + self.sortingInfo.fileInfo{td}.duration;
            end
            
            hold on;
            nSpikesUnit = histogram([Units.nSpikes], 100,'BinLimits', [0,10000]);
            barstrings = num2str(nSpikesUnit.Values);
             
            title('# Spikes per sec.');
            xlabel('# of spikes / sec');
            ylabel('# of units');
%            xl = xlim;
%            yl = ylim;
            maxText2 = sprintf('mean # of spikes');
            hold off;
            
            
            %%  Saving figures
            %  ############################################################
            
            if ~isempty(P.savingFolder)
                if exist(P.savingFolder)
                    disp('final output folder already exists')
                else
                    mkdir(P.savingFolder)
                    disp('final output folder created')
                end
                
                savefig(F.fh, fullfile(P.savingFolder, strcat(self.name,'-QC.fig')));
                saveas(F.fh, fullfile(P.savingFolder, strcat(self.name,'-QC.png')),'png');
                close all
            end            
            
        end
        
        % -----------------------------------------------------------------
        function FP = plotFootprints(self, varargin)
            P.Units = self.Units;
            P.scaling = [];
            P.wfsButtonDownFcn = [];
            P.maxNumberOfChannels = 10;
            P.fh = '';
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            Units = self.getUnits(P.Units);
            
            if numel(Units) == 1
                FP = Units.plotFootprint('IDs', Units.ID, ...
                    'scaling', P.scaling, ...
                    'wfsButtonDownFcn', P.wfsButtonDownFcn);
            else
                if isempty(self.MultiElectrode)
                    self.MultiElectrode = self.DataSource.MultiElectrode;
                end
                ME = self.MultiElectrode    
                
                fp = self.getFootprint(Units);
                FP = hdsort.plot.Waveforms2D(fp, ME.electrodePositions, ...
                    'maxNumberOfChannels', P.maxNumberOfChannels,...
                    'scaling', P.scaling, ...
                    'wfsButtonDownFcn', P.wfsButtonDownFcn, ...
                    'fh', P.fh);
            end
        end
        
        % -----------------------------------------------------------------
        function savedObject = saveobj(self)
            savedObject = self;
            % Nothing to do yet
        end
        
        % -----------------------------------------------------------------
        function savedObject = returnAsStruct(self)
            
            %% Turn all subfields either into structures or remove them:
            savedObject = struct(self);
            for ii = 1:numel(savedObject.Units)
                u_ = struct(savedObject.Units);
                u_ = rmfield(u_, 'parentPopulation');
                u_ = rmfield(u_, 'registeredAnalysis');
                u_ = rmfield(u_, 'plotManager');

                Units(ii) = u_;
                Units(ii).MultiElectrode = struct(Units(ii).MultiElectrode);
                Units(ii).MultiElectrode = rmfield(Units(ii).MultiElectrode, 'dataSource');
            end
            
            savedObject.Units = Units;
            savedObject.MultiElectrode = struct(savedObject.MultiElectrode);
            savedObject.MultiElectrode = rmfield(savedObject.MultiElectrode, 'dataSource');
            savedObject = rmfield(savedObject, 'DataSource')
        end
    end
    
    methods (Static)
        function savedObject = loadobj(savedObject)
            % Nothing to do yet
            savedObject.fileLocation = hdsort.hdsort.util.convertPathToOS(savedObject.fileLocation);
            savedObject.file = hdsort.hdsort.util.convertPathToOS(savedObject.file);
        end
    end
    
end