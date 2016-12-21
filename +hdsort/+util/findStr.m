function idx = findStr(stringCell, string)
% This function returns the location of a specific string in a cell of
% strings.
idx = find(~cellfun(@isempty, strfind(stringCell, string)));