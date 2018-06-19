function fileWrappers = checkFileWrappers(fileWrappers)
if ~iscell(fileWrappers)
    fileWrappers = {fileWrappers};
end
for i=1:length(fileWrappers)
    if ~isempty(fileWrappers{i})
        if isnumeric(fileWrappers{i})
            fileWrappers{i} = hdsort.file.DataMatrix(fileWrappers{i});
        else
            assert(isa(fileWrappers{i}, 'hdsort.file.FileWrapperInterface'), 'all datasources must implement the DataSourceInterface!');
        end
    end
end
