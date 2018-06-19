function M = matrixCreate(fname, h5path, dims, maxDims, h5type, chunkDims, deflation)
    if ~exist('chunkDims', 'var')
        chunkDims = []; deflation = [];
    end
    
    hdsort.filewrapper.hdf5.createVariableAndOrFile(fname, h5path, dims, maxDims, h5type, chunkDims, deflation);
    M = hdsort.filewrapper.hdf5.matrix(fname, h5path, false); 
end