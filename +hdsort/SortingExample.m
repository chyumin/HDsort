% 
% preprocessedFiles = {'some_file.h5', 'some_file2.h5'};
% DS = mysort.mea.CMOSMEA(preprocessedFiles);
% sortingPath = 'path-to-some-folder'
% sortingName = 'sortingNameXYZ';
% 
% %% Constructor:
% sorting = mysort.HDSorting.Sorting(DS, sortingPath, sortingName);
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
pd = pdefs()

rawFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'data', 'Trace_20160513_12_07_48.raw.h5')
configFile = fullfile( pd.mea1kRoot, 'rolandd', '160513', 'config', '11_30_06.mapping.nrk')
preprocessedFolder = fullfile(pd.mea1kIntermediate, 'rolandd', '161220');
sortingName = 'testGrid_161220'

%%
preprocessjob = grid.PreprocessJob(sortingName, preprocessedFolder, {rawFile}, 'configFile', configFile);
if ~preprocessjob.alreadyPreprocessed()
    preprocessjob.setTaskParameters();
    preprocessjob.createAutoSubmitToken();
    all_tasks_completed = preprocessjob.waitForTasksToFinish(10);
end

[preprocessedFiles, rawFiles, cmdFiles] = preprocessjob.getFileNames();
%%

DS = mysort.mea.CMOSMEA(preprocessedFiles);
sorting = mysort.HDSorting.Sorting(DS, preprocessedFolder, sortingName)
%%
sorting.startSorting('sortingMode', 'BSSE')
%%
sorting.startSorting('sortingMode', 'grid_for')
%%
sorting.postprocessGridSorting()
