function b = fileWrapperIsEmpty(fw)
    b = true;
    if isempty(fw)
        return
    end
    if iscell(fw)
        b = isempty(fw{1});
    end
        