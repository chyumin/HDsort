% -------------------------------------------------------------------------
%% New test scenario:
pd = hdsort.pathDefinitions()

if 0
    rawFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'data', 'Trace_20160513_12_07_48.raw.h5')
    configFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'config', '11_30_06.mapping.nrk')
    mainFolder = fullfile(pd.mea1kIntermediate, 'rolandd', '161220');
    sortingName = 'testGrid_161220'
    preprocessjob = hdsort.grid.PreprocessJob(sortingName, mainFolder, {rawFile}, 'configFile', configFile);
    
    if ~preprocessjob.alreadyPreprocessed()
        preprocessjob.setTaskParameters();
        preprocessjob.createAutoSubmitToken();
        all_tasks_completed = preprocessjob.waitForTasksToFinish(10);
    end
    
    [preprocessedFiles, rawFiles, cmdFiles] = preprocessjob.getFileNames();
    %%
else
    rawFile = fullfile('./Tutorial', 'new_file02.raw.h5')
    sortingName = 'tutorial_sorting';
    %mainFolder = fullfile(pd.mea1kIntermediate, 'rolandd', 'hd_tutorial');
    mainFolder = fullfile('./Tutorial')
    mkdir(mainFolder);
    RAW = hdsort.file.BELMEAFile(rawFile);
    preprocessedFiles = {RAW.preprocessFile(mainFolder)};
end

%%
DS = hdsort.file.CMOSMEA(preprocessedFiles);

%%
resultFile = fullfile(mainFolder, 'results.mat');

if exist(resultFile)
    load(resultFile)
else
    if 0
        sorting = hdsort.Sorting(DS, mainFolder, sortingName)
        %%
        [R, P] = sorting.startSorting('sortingMode', 'QSUB')
        %%
        [R, P] = sorting.startSorting('sortingMode', 'grid_for')
        
        %%
        [R, P] = sorting.startSorting('sortingMode', 'localHDSorting')
        
        %%
        [R, P] = sorting.postprocessGridSorting()
    else
        %% Alternate:
        [R, P] = hdsort.startHDSorting(DS, mainFolder, sortingName)
    end
    save(resultFile, 'R', 'P')
end

% -------------------------------------------------------------------------
%% Test hdsort.file.CortexLabFiles
folderSet = '/Volumes/hierlemann/intermediate_data/Mea1k/rolandd/cortexlab/set1';
datFile_ = dir([folderSet '/*.dat']);
datFile = fullfile(folderSet, datFile_.name);
% Read PRM file:
prmFile_ = dir([folderSet '/*.prm']);
prmFile = fullfile(folderSet, prmFile_.name);
%
M = hdsort.file.CortexLabFile(datFile, prmFile)
%%
preM = M.preprocessFile('/Users/rolandd/tmp')

% -------------------------------------------------------------------------
%% Test hdsort.file.CMOSMEA
folder = '/Volumes/hierlemann/intermediate_data/Mea1k/rolandd/160513/preprocessed';
fileList = {[folder '/Trace_20160513_11_32_04.h5'], [folder '/Trace_20160513_11_34_21.h5']}
mainFolder = folder;

%%
%DS = hdsort.file.CMOSMEA(fileList{1});
%covest_old = DS.getCovest();


%% Test single file hdsort.file.CMOSMEAFile
%cmosmeafile = hdsort.file.CMOSMEAFile(fileList{1});
%clist = cmosmeafile.getChannelList();
%covest = cmosmeafile.getCovest();

%% Test multifile hdsort.file.CMOSMEA
cmoslist = hdsort.file.CMOSMEA(fileList);
%covest_list = cmoslist.getCovest();

%%
cmoslist.restrictToChannels(1:2);
wfs1 = cmoslist.getWaveform(100, 20, 50);

cmoslist.restrictToChannels();
cmoslist.restrictToChannels(2:3);
wfs2 = cmoslist.getWaveform(100, 20, 50);



%%
DSList.restrictToChannels();
cmoslist.restrictToChannels();

DSList.restrictToChannels(1:10);
cmoslist.restrictToChannels(1:10);

%%
nstdold = DSList.noiseStd();
nstdnew = cmoslist.noiseStd();


%%
stold = DSList.detectSpikes('thr', 5);
stnew = cmoslist.detectSpikes('thr', 5);

%%
cutLeft = 20;
cutLength = 50;
channels = 1:10;

wfold = DSList.getWaveform(stold{1}, cutLeft, cutLength, channels);
wfnew = cmoslist.getWaveform(stold{1}, cutLeft, cutLength, channels);

%%
sortingName = 'asdf0123'
sortingPath = '/Volumes/hierlemann/intermediate_data/Mea1k/rolandd/170104'
sorting = hdsort.Sorting(cmoslist, sortingPath, sortingName)

%%
sorting.startSorting('sortingMode', 'localHDSorting')

%%
sorting.startSorting('sortingMode', 'grid_for')
sorting.postprocessGridSorting()
%%
sorting.startSorting('sortingMode', 'QSUB')

%%
dc = diag(covest.CCol);
dc_old = diag(covest_old.CCol);

figure; plot(dc); hold on; plot(dc_old+1000)


%%
M = randn(10000, 100);
MM = hdsort.file.DataMatrix(M);

%%
MM.restrictToChannels(1:2);
wfs1 = MM.getWaveform(100, 20, 50);
st1 = MM.detectSpikes();
x1 = MM.getData(1:10,:);

MM.restrictToChannels();
MM.restrictToChannels(2:3);
wfs2 = MM.getWaveform(100, 20, 50);
st2 = MM.detectSpikes();
x2 = MM.getData(1:10,:);

