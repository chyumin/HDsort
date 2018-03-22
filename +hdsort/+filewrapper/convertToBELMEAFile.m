function fileName = convertToBELMEAFile(fileName, dataMatrix, frameNumbers, ME, channelIdx)

delete(fileName)

if nargin < 4
    channelIdx = 1:size(dataMatrix, 2);
end

%% Set File as being in process
proc = hdsort.filewrapper.hdf5.createVariableAndOrFile(fileName, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
proc(1,1) = int32(1);

%% Set chipinformation:
hdf5write(fileName,'/chipinformation/software_version', '20170126', 'WriteMode', 'append')
hdf5write(fileName,'/chipinformation/chip_version', 'Mea1k', 'WriteMode', 'append')

%% Save mapping
nC_effective = numel(ME.electrodeNumbers);
assert(numel(channelIdx) == nC_effective, 'Effective number of channels must correspond to MultiElectrode')

mapping = struct();
for ii = 1:nC_effective
    mapping(ii).electrode = ME.electrodeNumbers(ii);
    mapping(ii).x = ME.electrodePositions(ii, 1);
    mapping(ii).y = ME.electrodePositions(ii, 2);
    mapping(ii).channel = ii-1;
end
hdf5write(fileName,'/ephys/mapping', mapping, 'WriteMode', 'append')

%% Save frame_rate
hdf5write(fileName,'/ephys/frame_rate', 20000, 'WriteMode', 'append')

%% Save framenumbers:
[nSamples, nC_original] = size(dataMatrix);
frameNumbers = frameNumbers(:);
assert( nSamples == size(frameNumbers, 1), 'Number of frames must be the same as the number of samples!')
frames = hdsort.filewrapper.hdf5.createVariableAndOrFile(fileName, '/ephys/frame_numbers', ...
    [nSamples 1], [nSamples 1], 'H5T_NATIVE_ULONG');
frames(:, 1) = frameNumbers;

%% Save signal:
dims = [nSamples nC_effective];
maxDims = [nSamples nC_effective];
h5Type = 'H5T_NATIVE_USHORT';
signal = hdsort.filewrapper.hdf5.createVariableAndOrFile(fileName, '/ephys/signal', ...
    dims, maxDims, h5Type); %, chunkDims, P.deflation);

chunkSize = 250000;
C = hdsort.util.Chunker(nSamples, 'chunkSize', chunkSize, 'progressDisplay', 'console');
while C.hasNextChunk()
    edges = C.getNextChunk();
    signal(edges(1):edges(2),:) = dataMatrix(edges(1):edges(2), channelIdx);
end

%%
clear signal
clear frames

% Set File as being done
proc(1,1) = int32(0);
clear proc
end