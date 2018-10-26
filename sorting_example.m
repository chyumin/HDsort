%% Add 'External' functions to path:
setup()

%% Download test file:
mainFolder = fullfile('.', 'Tutorial');
fileName = 'chirpsweep01.raw.h5';
rawFile = fullfile(mainFolder, fileName);
download_testfile(rawFile)
assert(exist(rawFile, 'file')>0, 'This tutorial requires a raw data file!')

%% Get a handle onto the raw file, which will then be used as input to the sorter:
sortingName = 'sorting_example';
RAW = hdsort.file.BELMEAFile(rawFile);
RAW.restrictToConnectedChannels(); % This line is very important to not sort empty electrodes!

%% Create the object that performs the sorting:
HDSorting = hdsort.Sorting(RAW, mainFolder, sortingName)

%% Preprocess:
% The preprocessor loads data from the file in chunks, filters it, and saves 
% the filtered data into a new hdf5 file that is standardized for all types
% of input data.
% If further performs a couple of operation for each local electrode group
% (LEG) such as spike detection, spike waveform cutting and noise
% estimation. The result is a folder named group000x for each LEG that
% contains the data necessary to perform the parallel parts of the sorting.
% This implementation minimizes the number of time-consuming file access
% operations and thus speeds up the parallel processes significantly, even
% allowing the parallel processes to be run on a small desktop computer or
% laptop.

chunkSize = 5e5; % This number depends a lot on the available RAM
HDSorting.preprocess('chunkSize', chunkSize, 'forceFileDeletionIfExists', true);

%% Sort each LEG independently:
HDSorting.sort('sortingMode', 'local_parfor'); % (default)
% Alternative sorting modes are:
% HDSorting.sort('sortingMode', 'local'); % for loop over each LEG
% HDSorting.sort('sortingMode', 'grid'); % requires a computer grid architecture

%% Combine the resutls of each LEG in the postprocessing step:
HDSorting.postprocess()

%% Export the results in an easy to read format:
[sortedPopulation, sortedPopulation_discarded] = HDSorting.createSortedPopulation(mainFolder)

%% When the sorting has already been run before, open the results from a file with:
sortedPopulation = hdsort.results.Population(HDSorting.files.results)

%% Plot some results right away:
sortedPopulation.plot();
sortedPopulation.Units(10).plot();


