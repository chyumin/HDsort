%%
if 0
    rawFile = fullfile('.', 'Tutorial', 'file01.raw.h5')
    RAW = hdsort.filewrapper.BELMEAFile(rawFile)
    me = RAW.getMultiElectrode()
    
    %%
    preprocessLocation = fullfile('.', 'Tutorial');
    preprocessedFile = RAW.preprocessFile(preprocessLocation)
    preprocessedFiles = {preprocessedFile}
    
    %%
    sortingName = 'hd_sort_tutorial01';
    rawFiles = {rawFile};
    preprocessjob = hdsort.grid.PreprocessJob(sortingName, preprocessLocation, rawFiles, 'gridType', 'QSUB')
    preprocessjob.alreadyPreprocessed()
    if ~preprocessjob.alreadyPreprocessed()
        preprocessjob.setTaskParameters();
        preprocessjob.createAutoSubmitToken();
        all_tasks_completed = preprocessjob.waitForTasksToFinish(10);
    end
    [preprocessedFiles, rawFiles, cmdFiles] = preprocessjob.getFileNames();
    
end

baseFolder = './Tutorial';

%% 1. Reshape the original recordings such that they contain two interchangable recording areas:
original_rawFile = fullfile('.', 'Tutorial', 'file02.raw.h5')
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

new_rawFile = fullfile('.', 'Tutorial', 'new_file02.raw.h5')
if ~exist(new_rawFile)
    new_rawFile = hdsort.filewrapper.convertToBELMEAFile(new_rawFile, dataMatrix, frameNumbers, newME, channelIdx)
end

%% Preprocess this file:
newRAW = hdsort.filewrapper.BELMEAFile(new_rawFile)
preprocessLocation = fullfile('.', 'Tutorial');
preprocessedFile = newRAW.preprocessFile(preprocessLocation)
preprocessedFiles = {preprocessedFile}

%% Sort original recordings in order to get neuron candidates:
PRE = hdsort.filewrapper.CMOSMEA(preprocessedFiles);
sortingLocation =  fullfile('.', 'Tutorial');
sorting = hdsort.Sorting(PRE, sortingLocation, sortingName)
[R, P] = sorting.startSorting('sortingMode', 'localHDSorting')

%% Create artificialUnits:
rawFileNames = {new_rawFile};

%%
hdbenchmarking.generate.artificialUnits2






