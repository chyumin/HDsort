function Y = readBinDtype(fname, DTYPE)
 
    fh = fopen(fname,'r');
    try
        [filename, permission, machineformat, encoding] = fopen(fh);
        i = 1;
        idx = [];
        for k=1:size(DTYPE,1)
            if strcmp('skip', DTYPE{k,2})
                continue
            end
            idx(end+1) = k;

            offset = sum(cell2mat(DTYPE(1:k-1,3)));
            fseek(fh, offset,'bof');
            skip = sum(cell2mat(DTYPE(setdiff(1:size(DTYPE,1),k),3)));
                  
            Y.(DTYPE{k,1}) = fread(fh, inf, DTYPE{k,2}, skip);
            i = i+1;
        end
        fclose(fh);
    catch
        fclose(fh);
        rethrow(lasterror);
    end
    %Y = cell2struct(X,DTYPE(idx,1),2);