function fileWrappers = checkFileWrappers(fileWrappers)
if ~iscell(fileWrappers)
    fileWrappers = {fileWrappers};
end
for i=1:length(fileWrappers)
    if ~isempty(fileWrappers{i})
        if isnumeric(fileWrappers{i})
            fileWrappers{i} = mysort.ds.Matrix(fileWrappers{i});
        else
            assert(isa(fileWrappers{i}, 'hdsort.filewrapper.FileWrapperInterface'), 'all datasources must implement the DataSourceInterface!');
        end
    end
end
