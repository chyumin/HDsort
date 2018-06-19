%%
baseFolder = './Tutorial';

% -------------------------------------------------------------------------
%% 1. Reshape the original recordings such that they contain two interchangable recording areas:
original_rawFile = fullfile('.', 'Tutorial', 'example_data', 'file01.raw.h5')
originalRAW = hdsort.filewrapper.BELMEAFile(original_rawFile)
originalME = originalRAW.getMultiElectrode()

% Manually select the index of one electrode in each block that correspond to the
% location in the other block. They will be the anchor to overlap the
% blocks:
originalME.plotConfig();
el1 = 515;
el2 = 589;

%% Generate new ME:
[newME, swapElectrodePairs, blockIdx] = hdbenchmarking.generate.overlappingBlocks(originalME, el1, el2, 1);
[~, name_] = fileparts(original_rawFile); [~, fileBaseName] = fileparts(name_);
newMultiElectrodeFile = fullfile(baseFolder, [fileBaseName '_newMultiElectrode.mat']);
save(newMultiElectrodeFile, 'newME', 'swapElectrodePairs', 'blockIdx', 'fileBaseName');

%% Generate a new file:
% To create a file in the BELMEA file format, you need to create a matrix object
% containing the raw data (either using the hdf5.matrix wrapper or for any
% matlab matrix the filewrapper.DataMatrix). You further need the frameNumbers,
% the new MultiElectrode and an index that matches the channels in the
% original dataMatrix with the new MultiElectrode.
dataMatrix = hdsort.filewrapper.hdf5.matrix( originalRAW.fileName, '/ephys/signal')
frameNumbers = originalRAW.getFrameNumbers();

channelIdx_ = false(size(dataMatrix, 2), 1);
for ii = 1:numel(newME.electrodeNumbers)
    c = find(newME.electrodeNumbers(ii) == originalME.electrodeNumbers);
    channelIdx_(c) = true;
end
channelIdx = find(channelIdx_);

new_rawFile = fullfile('.', 'Tutorial', 'new_file01.raw.h5')
if ~exist(new_rawFile)
    new_rawFile = hdsort.filewrapper.convertToBELMEAFile(new_rawFile, dataMatrix, frameNumbers, newME, channelIdx)
end

% -------------------------------------------------------------------------
%% 2. Sort original recordings in order to get neuron candidates:
sortingName = 'hd_sort_tutorial01'
sortingLocation =  fullfile('.', 'Tutorial');

newRAW = hdsort.filewrapper.BELMEAFile(new_rawFile)

sorting = hdsort.Sorting(newRAW, sortingLocation, sortingName)
sorting.preprocess('forceFileDeletionIfExists', 1)
sorting.sort('sortingMode', 'local')
sorting.postprocess()

% Create a SpikeSortingResult
SpikeSortingResult = sorting.createSpikeSortingResult(sortingLocation)
nSpikes = SpikeSortingResult.getSpikeCounts();

% -------------------------------------------------------------------------
%% 3. Generate artificial data (example)
baseFolder = './Tutorial';

rawFileName = fullfile(baseFolder, 'new_file01.raw.h5');
preFileNames = sorting.files.preprocessed;

% Select the electrodes that demark the corners of the hidens blocks by
% hand:
el1 = 381;
el2 = 382;

resultFile = fullfile(baseFolder, 'results.mat');
load(resultFile)

datasetName = 'amplitude_sweep01';
gdf = double(R.gdf_merged);
estimated_footprints =  R.T_merged;
RAW_list = {hdsort.filewrapper.BELMEAFile(rawFileName)};
PRE = hdsort.filewrapper.CMOSMEA(preFileNames);
[~, swapElectrodePairs, blockIdx] = hdbenchmarking.generate.overlappingBlocks(PRE.MultiElectrode, el1, el2, true);

[artificialFileList, artificialUnitFile] = hdbenchmarking.generate.artificialUnits(datasetName, ...
    baseFolder, gdf, estimated_footprints, RAW_list, PRE, swapElectrodePairs, blockIdx)

% -------------------------------------------------------------------------
%% 4. Sort the artificial dataset:
artificialPRE = hdsort.filewrapper.CMOSMEA(artificialFileList);
artificialSortingName = 'hdbenchmarking_tutorial01';

sortingLocation =  fullfile('.', 'Tutorial');

sorting_artificial = hdsort.Sorting(artificialPRE, sortingLocation, artificialSortingName)
sorting_artificial.preprocess('forceFileDeletionIfExists', 1, 'prefilter', 0)
sorting_artificial.sort('sortingMode', 'local')
sorting_artificial.postprocess()

SpikeSortingResult_artificial = sorting_artificial.createSpikeSortingResult(sortingLocation)

% -------------------------------------------------------------------------
%% 5. Analyze the results by comparing them to the ground truth:
analysisLocation =  fullfile('.', 'Tutorial');
artificialUnits = load(artificialUnitFile)
[sortingEvaluation, sortingEvaluationFile] = hdbenchmarking.evaluate.sorting(...
    artificialUnits, SpikeSortingResult_artificial, analysisLocation);

% -------------------------------------------------------------------------
%% 6. Plot the results:
figure; subplot(3,1,1)
scatter(sortingEvaluation.gt.AMP, sortingEvaluation.matched.TPR)
ylabel('sensitivity (TPR) [%]')

subplot(3,1,2)
scatter(sortingEvaluation.gt.AMP, sortingEvaluation.matched.PPV)
ylabel('precision (PPV) [%]')

subplot(3,1,3)
scatter(sortingEvaluation.gt.AMP, sortingEvaluation.matched.ERR)
ylabel('error rate [s^(-1)]')




