function [PopulationResult, PopulationResult_discarded, P] = ...
    createPopulation(sortingName, results, rawDSList, noiseStd, varargin)

P.analysisFolder = '';
P = hdsort.util.parseInputs(P, varargin, 'error');

disp('Creating Population:')

%% Turn results into double format:
gdf_merged = double(results.gdf_merged);
footprints_merged = double(results.T_merged);
gdf_discarded = double(results.gdf_discarded);
footprints_discarded = double(results.T_discarded);

%% Break gdfs apart:
L = cellfun(@(x) size(x, 1), rawDSList);
trialTimes = [0 cumsum(L)];
assert(max(gdf_merged(:,2)) <= trialTimes(end), 'Spike Times out of Range !!');
[mgdf_merged, nS_merged] = hdsort.spiketrain.gdf2multiSessionGdf(gdf_merged, trialTimes);
[mgdf_discarded, nS_discarded] = hdsort.spiketrain.gdf2multiSessionGdf(gdf_discarded, trialTimes);

if any(round(mgdf_merged(:,2)) ~= mgdf_merged(:,2))
    warning('Spiketimes had to be rounded!!!')
    mgdf_merged(:,2) = round(mgdf_merged(:,2));
    
    if ~isempty(mgdf_discarded)
        mgdf_discarded(:,2) = round(mgdf_discarded(:,2));
    end
end
startTimes = [];

%% Get Multielectrode
nFiles = numel(rawDSList);
MultiElectrode = rawDSList{1}.MultiElectrode;
for fi = 2:nFiles
    assert(MultiElectrode == rawDSList{fi}.MultiElectrode, 'The raw files do not share the same Multielectrode!')
end

%% Go through each raw file and create some statistics:
fileInfo = {}; missingFrames = {};

startSpikeIdx = 1; lastSpikeIdx = 0;
startSpikeIdx_discarded = 1; lastSpikeIdx_discarded = 0;
for fi = 1:nFiles
    disp(['Align the sorting output with the original frame numbers. ' num2str(fi) ' out of ' num2str(nFiles)] );
    
    Fs = rawDSList{fi}.samplesPerSecond;
    FN = rawDSList{fi}.getFrameNumbers();
    
    fileInfo{fi}.missingFrames = rawDSList{fi}.getMissingFrameNumbers();
    fileInfo{fi}.duration = (FN(end)-FN(1))/Fs;
    fileInfo{fi}.nSamples = numel(FN);
    fileInfo{fi}.firstFrame = FN(1);
    fileInfo{fi}.lastFrame = FN(end);
    startTimes(fi,:) = [FN(1), FN(end)];
    
    %% Convert gdf-spiketimes relative to file beginning to absolute spiketimes:
    lastSpikeIdx = lastSpikeIdx + nS_merged(fi);
    gdf_merged(startSpikeIdx:lastSpikeIdx, 2) = FN( mgdf_merged( startSpikeIdx:lastSpikeIdx , 2) );
    startSpikeIdx = lastSpikeIdx + 1;
    
    if ~isempty(nS_discarded)
        lastSpikeIdx_discarded = lastSpikeIdx_discarded + nS_discarded(fi);
        gdf_discarded(startSpikeIdx_discarded:lastSpikeIdx_discarded, 2) =  FN( mgdf_discarded( startSpikeIdx_discarded:lastSpikeIdx_discarded , 2));
        startSpikeIdx_discarded = lastSpikeIdx_discarded + 1;
    end
end
clear 'FN'

%% Create Population:
if ~isfield(results, 'cutLeft')
    warning('This should not happen anymore!')
    try
        results.cutLeft = results.G(1).templates.cutleft;
        for G = results.G
            assert( ~any( G.templates.cutleft - results.cutLeft), 'Some templates have different cutLefts!')
        end
    catch
        results.cutLeft = [];
        warning('Missing cutLeft in sorting results!')
    end
end


disp('Create hdsort.results.Population...');
sortingInfo.fileInfo = fileInfo;
sortingInfo.startTimes = startTimes;
sortingInfo.summary = results.summary;
sortingInfo.sortingParameters = results.sortingParameters;

PopulationResult = hdsort.results.Population(...
    'gdf', gdf_merged, ...
    'footprints', footprints_merged, ...
    'name', sortingName, ...
    'noiseStd', noiseStd, ...
    'fileLocation', P.analysisFolder, ...
    'MultiElectrode', MultiElectrode, ...
    'cutLeft', results.cutLeft, ...
    'sortingInfo', sortingInfo);

if ~isempty(gdf_discarded)
    PopulationResult_discarded = hdsort.results.Population(...
        'gdf', gdf_discarded, ...
        'footprints', footprints_discarded, ...
        'name', [sortingName '_discarded'], ...
        'noiseStd', noiseStd, ...
        'fileLocation', P.analysisFolder, ...
        'MultiElectrode', MultiElectrode, ...
        'cutLeft', results.cutLeft, ...
        'sortingInfo', sortingInfo);
else
    PopulationResult_discarded = [];
end

end
