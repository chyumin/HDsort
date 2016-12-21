function str = dlmstring(M, delim, format)
    if ~iscell(M)
        M = num2cell(M);
    end
    
    if isempty(M)
        str = [];
        return
    end
    
    % creates a delimited string
    if ~exist('delim', 'var')
        delim = ',';
    end
    if ~exist('format', 'var')
        format = '%d';
    end

    str = sprintf(format, M{1});
    for i=2:length(M)
        str = [str delim sprintf(format, M{i})];
    end