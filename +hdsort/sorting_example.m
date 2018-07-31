% -------------------------------------------------------------------------
%% New test scenario:
setup()
pd = hdsort.pathDefinitions()


%% Download test file:
mainFolder = fullfile('.', 'Tutorial');
    
fileName = 'chirpsweep01.raw.h5';
rawFile = fullfile(mainFolder, fileName);

if ~exist(rawFile)
    
    userInput = input(['Do you want to download test file to ' mainFolder ' (~2GB)? Y/N:'],'s');
    assert(strcmp(userInput, 'Y'), 'User denied')

    mkdir(mainFolder);
    url = 'https://polybox.ethz.ch/index.php/s/4JqS6b0q2VGZKs6/download'
    websave(rawFile,url)
end
assert(exist(rawFile)>0, 'This tutorial requires a raw data file!')

%%
sortingName = 'sorting_example';
RAW = hdsort.file.BELMEAFile(rawFile);
RAW.restrictToConnectedChannels(); % This line is very important to not sort empty electrodes!

HDSorting = hdsort.Sorting(RAW, mainFolder, sortingName)

%%
chunkSize = 5e5; % This number depends a lot on the available RAM 
HDSorting.preprocess('chunkSize', chunkSize, 'forceFileDeletionIfExists', true);
%%
HDSorting.sort()
%HDSorting.sort('sortingMode', 'local')

%%
HDSorting.postprocess()

%%
[sortedPopulation, sortedPopulation_discarded] = HDSorting.createSortedPopulation(mainFolder)


