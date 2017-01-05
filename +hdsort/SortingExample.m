% 
% preprocessedFiles = {'some_file.h5', 'some_file2.h5'};
% DS = mysortx.mea.CMOSMEA(preprocessedFiles);
% sortingPath = 'path-to-some-folder'
% sortingName = 'sortingNameXYZ';
% 
% %% Constructor:
% sorting = mysortx.HDSorting.Sorting(DS, sortingPath, sortingName);
% 
% %% New Sorting:
% sorting.startSorting('runMode', 'grid', 'parfor',true);
% 
% %% Go back to a sorting that had been started on the grid before:
% sorting.reuptakeSorting();
% 
% %% Postprocessing:
% sorting.postprocessGridSorting();

% -------------------------------------------------------------------------
%% New test scenario:
pd = hdsort.pathDefinitions()

rawFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'data', 'Trace_20160513_12_07_48.raw.h5')
configFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'config', '11_30_06.mapping.nrk')
preprocessedFolder = fullfile(pd.mea1kIntermediate, 'rolandd', '161220');
sortingName = 'testGrid_161220'

%%
preprocessjob = hdsort.grid.PreprocessJob(sortingName, preprocessedFolder, {rawFile}, 'configFile', configFile);
if ~preprocessjob.alreadyPreprocessed()
    preprocessjob.setTaskParameters();
    preprocessjob.createAutoSubmitToken();
    all_tasks_completed = preprocessjob.waitForTasksToFinish(10);
end

[preprocessedFiles, rawFiles, cmdFiles] = preprocessjob.getFileNames();
%%

DS = hdsort.filewrapper.CMOSMEA(preprocessedFiles);
sorting = hdsort.Sorting(DS, preprocessedFolder, sortingName)
%%
sorting.startSorting('sortingMode', 'QSUB')
%%
sorting.startSorting('sortingMode', 'grid_for')
%%
sorting.postprocessGridSorting()


%% Alternate:
[gdf_merged, T_merged, localSorting, localSortingID, sessionLengths] = hdsort.startHDSorting(preprocessedFiles, preprocessedFolder, 'testGrid_startHDSorting')


% -------------------------------------------------------------------------
%% Test hdsort.filewrapper.CMOSMEA
folder = '/Volumes/hierlemann/intermediate_data/Mea1k/rolandd/160513/preprocessed';
fileList = {[folder '/Trace_20160513_11_32_04.h5'], [folder '/Trace_20160513_11_34_21.h5']}
preprocessedFolder = folder;

%%
%DS = hdsort.filewrapper.CMOSMEA(fileList{1});
%covest_old = DS.getCovest();

%%
DSList = mysortx.mea.CMOSMEA(fileList);
%covest_listold = DSList.getCovest();

%% Test single file hdsort.filewrapper.CMOSMEAFile
%cmosmeafile = hdsort.filewrapper.CMOSMEAFile(fileList{1});
%clist = cmosmeafile.getChannelList();
%covest = cmosmeafile.getCovest();

%% Test multifile hdsort.filewrapper.CMOSMEA
cmoslist = hdsort.filewrapper.CMOSMEA(fileList);
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
MM = hdsort.filewrapper.DataMatrix(M);

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

