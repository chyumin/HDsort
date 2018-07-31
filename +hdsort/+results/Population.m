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
        filePath
        file
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Population(sortingName, gdf_or_st, footprints, ...
             MultiElectrode, noiseStd, varargin)
            
            P.unitIDs    = [];
            P.DataSource = [];
            P.filePath = '';
            P.sortingInfo = [];
            P.cutLeft = [];
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            self.name = sortingName;
            self.noiseStd = noiseStd(:);
            self.MultiElectrode = MultiElectrode;
            
            self.filePath = P.filePath;
            self.sortingInfo = P.sortingInfo;
            
            %% Set datasource if it is provided as argument:
            if ~isempty(P.DataSource)
                self.DataSource = P.DataSource;
            else
                self.DataSource = [];
            end
            
            %% Initialize gdf, spikeTrains, unitIDs:
            if iscell(gdf_or_st)
                assert(length(P.unitIDs) == length(gdf_or_st), 'If spike trains are used to construct sorting, you MUST also give the unit IDs!');
                unitIDs = P.unitIDs;
                spikeTrains = gdf_or_st;
                gdf = hdsort.spiketrain.toGdf(spikeTrains, unitIDs);
            else
                assert(isempty(P.unitIDs), 'If input is a gdf, do not provide unitIDs!')
                s = size(gdf_or_st,2);
                assert(s >= 2 , 'Must have 2 or more columns!');
                gdf = gdf_or_st;
                [spikeTrains, unitIDs] = hdsort.spiketrain.fromGdf(gdf_or_st);
            end
            
            %% Initialize amplitudes if possible:
            if size(gdf_or_st,2) >= 2
                relAmpsV = gdf(:,3)./self.noiseStd(gdf(:,4));
            else
                relAmpsV = gdf(:, 1)*0 + 1;
            end
            
            assert(size(footprints,3) == length(unitIDs), 'There must be one footprint per Unit!');
            
            self.Units = hdsort.results.Unit.empty(0, length(unitIDs));
            for ui=1:length(unitIDs)
                IDs = unitIDs(ui);
                spikeTrain = spikeTrains{ui};
                
                if ~isempty(footprints)
                    footprint = footprints(:,:,ui);
                else
                    footprint = [];
                end
                
                self.Units(ui) = hdsort.results.Unit('ID', IDs, 'spikeTrain', spikeTrain, ...
                    'parentPopulation', self, ...
                    'footprint', footprint, 'cutLeft', P.cutLeft, ...
                    'spikeAmplitudes', relAmpsV(gdf(:, 1) == IDs), ...
                    'filePath', self.filePath);
            end
            
            self.getSpikeCounts();
        
        end % End of constructor
        
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
                folderName = fullfile(self.filePath, 'plots');
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
            savedObject.filePath = hdsort.util.convertPathToOS(savedObject.filePath);
            savedObject.file = hdsort.util.convertPathToOS(savedObject.file);
        end
    end
    
end