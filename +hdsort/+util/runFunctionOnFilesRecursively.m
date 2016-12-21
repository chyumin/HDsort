function runFunctionOnFilesRecursively(dpath, filefilter, fhandle)
    % recursively walks through dpath and runs fhandle on every file
    % defined by filefilter
    % filefilter has to be usable by dir, e.g. '2012_*.mat';
    
    fprintf('Processing directory: %s\n', dpath);
    allfnames = dir(dpath);
    % Run recursively 
    for i=1:length(allfnames)
        myFile = allfnames(i);
        if myFile.name(1) == '.'
            continue
        end
           
        if myFile.isdir
            hdsort.util.runFunctionOnFilesRecursively(fullfile(dpath, myFile.name), filefilter, fhandle);
        end
    end

    % Now run on files in this folder
    fnames = dir(fullfile(dpath, filefilter));
    for i=1:length(fnames)
        myFile = fnames(i);
        if myFile.name(1) == '.'
            continue
        end
           
        fhandle(fullfile(dpath, myFile.name));
    end        
end