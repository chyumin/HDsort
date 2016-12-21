
preprocessedFiles = {'some_file.h5', 'some_file2.h5'};
DS = mysort.mea.CMOSMEA(preprocessedFiles);
sortingPath = 'path-to-some-folder'
sortingName = 'sortingNameXYZ';

%% Constructor:
sorting = mysort.HDSorting.Sorting(DS, sortingPath, sortingName);

%% New Sorting:
sorting.startSorting('runMode', 'grid', 'parfor',true);

%% Go back to a sorting that had been started on the grid before:
sorting.reuptakeSorting();

%% Postprocessing:
sorting.postprocessGridSorting();


