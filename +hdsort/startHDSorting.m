function [R, P, P2] = startHDSorting(DSFull, outPath, runName, varargin)
% Spike Sorting of a High Density Configuration Array Recording
% This script can also be used to sort non-high density configurations
% but this might lead to unforseen problems.
%
% Input:
%   DSFull - Datasource Interface Instance of the Data. See
%            ana.startHDSortingTest.m for an example how it is used.
%   outPath - output folder where the sorting info should be stored
%   runName - (optional) Name of the run of this Sorting. If not
%             provided it defaults to "sortingRun".
useParForIfPossible = false;

% Sorting Parameters:
P = struct();
P.spikeDetection.method = '-';
P.spikeDetection.thr = 4.2;
P.artefactDetection.use = 0;
P.botm.run = 0;
P.spikeCutting.maxSpikes = 200000000000; % Set this to basically inf

P.noiseEstimation.minDistFromSpikes = 80;

P.spikeAlignment.initAlignment = '-';
P.spikeAlignment.maxSpikes = 60000;     % so many spikes will be clustered
P.clustering.maxSpikes = P.spikeAlignment.maxSpikes;  % dont align spikes you dont cluster...
P.clustering.meanShiftBandWidth = sqrt(1.8*6);
P.mergeTemplates.merge = 1;
P.mergeTemplates.upsampleFactor = 3;
P.mergeTemplates.atCorrelation = .93; % DONT SET THIS TOO LOW! USE OTHER ELECTRODES ON FULL FOOTPRINT TO MERGE
P.mergeTemplates.ifMaxRelDistSmallerPercent = 30;

P = hdsort.util.parseInputs(P, varargin, 'error');

name = 'dirtySort';
outPath = fullfile(outPath, runName);
if ~exist(outPath, 'file')
    mkdir(outPath);
end
if nargin < 3
    runName = 'sortingRun';
end

cutleft = 10;
Tf = 45;

MES = DSFull.MultiElectrode.toStruct(); % THIS LINE IS FOR SAVING !!
electrodePositions = DSFull.MultiElectrode.electrodePositions;
electrodeNumbers   = DSFull.MultiElectrode.electrodeNumbers;
% Make groupings of electrodes to be sorted independently
[groupsidx nGroupsPerElectrode] = hdsort.leg.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2));
% replace electrode indices with electrode numbers
groups = {};
for ii=1:length(groupsidx)
    groups{ii} = electrodeNumbers(groupsidx{ii});
end
groupFile = fullfile(fullfile(outPath, 'groupFile.mat'));
save(groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');

% Prepare data
useParFor = false;
if useParForIfPossible
    pPool = gcp;
    if isempty(pPool) pPool = parpool(); end
    if pPool.NumWorkers > 0 useParFor = true; end
end

if useParFor
    parfor ii=1:length(groupsidx)
        % Build Mea that concatenates multiple ntk files (WARNING NEED TO HAVE THE SAME CONFIGURATIONS!
        DScopy = hdsort.filewrapper.CMOSMEA(h5FileList, 'useFilter', 0, 'name', 'PREFILT');
        DScopy.restrictToChannels(groupsidx{ii});
        [S P_] = hdsort.sortScript(DScopy, fullfile(fullfile(outPath, 'groups'), ['group' sprintf('%04d', ii)]), runName, P);
        % RELEASE CHANNEL RESTRICTIONS FOR TEMPLATE ESTIMATION
        DScopy.restrictToChannels();
        hdsort.leg.startHDSortingTemplateEstimation(fullfile(outPath, 'groups'), ['group' sprintf('%04d', ii)], runName, Tf, cutleft, DScopy, groupsidx{ii}, MES);
    end
else
    for ii=1:length(groupsidx)
        DSFull.restrictToChannels(groupsidx{ii});
        [S P_] = hdsort.sortScript(DSFull, fullfile(fullfile(outPath, 'groups'), ['group' sprintf('%04d', ii)]), runName, P);
        % RELEASE CHANNEL RESTRICTIONS FOR TEMPLATE ESTIMATION
        DSFull.restrictToChannels();
        hdsort.leg.startHDSortingTemplateEstimation(fullfile(outPath, 'groups'), ['group' sprintf('%04d', ii)], runName, Tf, cutleft, DSFull, groupsidx{ii}, MES);
    end
end

clear groupDS S;

% Re-Estimate templates on all electrodes of that config after sorting
disp('DONE WITH ALL GROUPS')

%% Process all local sortings into a final sorting
groupFile = fullfile(fullfile(outPath, 'groupFile.mat'));
load(groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');

disp('Postprocessing...');
[R, P2] = hdsort.leg.processLocalSortings(outPath, runName, groups, groupsidx);
units = unique(R.gdf_merged(:,1));
nU = length(units)
assert(length(R.localSorting) == nU, 'must be identical');
assert(length(R.localSortingID) == nU, 'must be identical');
assert(size(R.T_merged,3) == nU, 'must be identical');
save(fullfile(outPath, [runName '_results.mat']), '-struct', 'R', '-v7.3');

