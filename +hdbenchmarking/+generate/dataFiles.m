function [newFileList, P] = generateArtificialData(sourceFileList, outputLocation, artificialUnitFile, varargin)
% The purpose of this script is to create a new file in the 'preprocessed'
% style taking a DataSource object as input and then adding neuron
% footPrints at given timepoints and with given amplitude.
P.chunkLen = [];
P.sessionName = '/Sessions/Session0';
P.save_as_binary = true;
P.debug = false;
P.deleteOldFilesIfNecessary = false;
P = util.parseInputs(P, varargin);

%% Check Inputs:
artificialUnits = load(artificialUnitFile);
auFileChecksum = Simulink.getFileChecksum(artificialUnitFile);

nFiles = numel(sourceFileList);
[nTf, nC, nU, nJ] = size(artificialUnits.fp_jittered);
[nU_, nFiles_] = size(artificialUnits.spikeTrains);
cutLeft = artificialUnits.cutLeft;

assert(nFiles == nFiles_, 'Not enough spiketrains!')
assert(nU == nU_, 'Not enough spiketrains!')
assert(cutLeft < nTf, 'Cut left must be shorter than the entire footprint!')

%%
newFileList = {};
for fi = 1:nFiles
    
    %% Set new file up:
    sourceFile = sourceFileList{fi};
    newFileList{fi} = fullfile(outputLocation, ['file' num2str(fi, '%02i') '_' artificialUnits.datasetName '.h5']);
    newFileName = newFileList{fi};
    
    [pathstr,name,ext] = fileparts(newFileName);
    binFile = fullfile(pathstr, [name, '.dat']);
    
    %% Delete file if necessary:
    if P.debug
        delete(newFileName)
        [pathstr,name,ext] = fileparts(newFileName);
        binFile = fullfile(pathstr, [name, '.dat']);
        delete(binFile);
    end
    
    %% Check whether the checksum of the artificial correspond to this file, if it already exists:
    if exist(newFileName)
        try
            auCS = h5read(newFileName, '/artificialUnits/checksum');
            auCS = auCS{1}(1:end-1);
        catch
            auCS = '';
        end
        
        if ~P.deleteOldFilesIfNecessary
            assert(strcmp(auCS, auFileChecksum), 'The existing file was not created with the same artificial units!')        
            disp('File has already been created with these settings')
            continue;
        else
            if ~strcmp(auCS, auFileChecksum)
                warning('Old file with different artificial units deleted!')
                delete(newFileName)
                delete(binFile)
            else
                continue;
            end
        end
    end
    
    %%
    DS = mysort.mea.CMOSMEA(sourceFile);
    ME = DS.MultiElectrode;
    [nSamples, nC_] = size(DS);
    assert(nC_ == nC, 'Channel-numbers must correspond!');
    
    %% Set File as being in process
    assert(~exist(newFileName, 'file'), ['Output File does already exist! ' newFileName]);
    
    proc = mysort.h5.createVariableAndOrFile(newFileName, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
    proc(1,1) = int32(1);
    
    %% Save checksum of artificial units:
    hdf5write(newFileName,'/artificialUnits/checksum', auFileChecksum, 'WriteMode', 'append')
    
    %% Set filter info:
    pref = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/prefiltered'], [1 1], [1 1], 'H5T_NATIVE_INT');
    pref(1,1) = h5read(sourceFile, [P.sessionName '/filter/prefiltered']);
    clear pref
    high = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/highpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
    high(1,1) = h5read(sourceFile, [P.sessionName '/filter/highpass']);
    clear high
    low = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/lowpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
    low(1,1) = h5read(sourceFile, [P.sessionName '/filter/lowpass']);
    clear low
    down = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/downsamplefactor'], [1 1], [1 1], 'H5T_NATIVE_INT');
    down(1,1) = h5read(sourceFile, [P.sessionName '/filter/downsamplefactor']);
    clear down
    ord = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/order'], [1 1], [1 1], 'H5T_NATIVE_INT');
    ord(1,1) =  h5read(sourceFile, [P.sessionName '/filter/order']);
    clear type
    ord = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/type'], [1 20], [1 20], 'H5T_C_S1');
    filterType = h5read(sourceFile, [P.sessionName '/filter/type']); filterType = [filterType{:}];
    ord(1,1:length(filterType)) = filterType;
    clear ord
    gd = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/filter/gainmultiplier'], [1 1], [1 1], 'H5T_NATIVE_INT');
    gainmultiplier = double(h5read(sourceFile, [P.sessionName '/filter/gainmultiplier']));
    gd(1,1) = 1; %gainmultiplier;
    clear gd
    
    % SAVE THE SOURCEFILES
    ffile = h5read(sourceFile, [P.sessionName '/source_files/raw_h5']); ffile = [ffile{:}];
    ord = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/source_files/raw_h5'], [1 length(ffile)], [1 length(ffile)], 'H5T_C_S1');
    ord(1,1:length(ffile)) = ffile;
    clear ord
    %mfile = h5read(sourceFile, [P.sessionName '/source_files/mapping_file']); mfile = [mfile{:}];
    %ord = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/source_files/mapping_file'], [1 length(mfile)], [1 length(mfile)], 'H5T_C_S1');
    %ord(1,1:length(mfile)) = mfile;
    %clear ord
    
    % CHIP ID
    chipid = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/chipid'], [1 1], [1 1], 'H5T_NATIVE_INT');
    chipid(1,1) =  h5read(sourceFile, [P.sessionName '/chipid']);
    clear chipid
    
    % ADC range and resolution
    gain = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/adc_resolution'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
    gain(1,1) = 1;%h5read(sourceFile, [P.sessionName '/adc_resolution']);
    gain = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/adc_range'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
    gain(1,1) = 1;%h5read(sourceFile, [P.sessionName '/adc_range']);
    clear gain;
    
    % VERSION
    version = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/version'], [1 1], [1 1], 'H5T_NATIVE_INT');
    version(1,1) =  h5read(sourceFile, [P.sessionName '/version']);
    clear version
    
    % GAIN
    gain = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/gain'], [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
    gain(1,1:4) = h5read(sourceFile, [P.sessionName '/gain']);
    clear gain
    
    % SR
    sr_ = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/sr'], [1 1], [1 1], 'H5T_NATIVE_INT');
    sr_(1,1) = h5read(sourceFile, [P.sessionName '/sr']); %int32(DS.getSampleRate);
    clear sr_
    
    % CHANNEL LIST
    names = {'channel_nr', 'connected', 'x', 'y', 'idx', 'dummy', 'damaged'};
    type_id = H5T.create('H5T_COMPOUND',length(names)*32);
    for i=1:length(names)
        H5T.insert(type_id, names{i}, (i-1)*32, 'H5T_NATIVE_INT');
    end
    cl = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/channel_list'], nC, nC, type_id);
    
    x = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/channel_nr'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
    x(1,1:nC) = ME.electrodeNumbers;
    clear x
    x = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/channel_posx'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
    x(1,1:nC) = ME.electrodePositions(:,1)';
    clear x
    x = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/channel_posy'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
    x(1,1:nC) = ME.electrodePositions(:,2)';
    clear x
    x = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/channel_connected'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
    x(1,1:nC) = ones(1,nC);
    clear x
    
    
    if isempty(P.chunkLen) || (~isempty(P.chunkLen) && P.chunkLen == 0)
        % This disables chunking and deflation
        chunkDims = [];
    else
        if P.chunkIndividualChannels
            chunkDims = [P.chunkLen 1];
        else
            chunkDims = [P.chunkLen nC];
        end
    end
    lastSample = 0;
    
    dims = [nSamples nC];
    maxDims = [nSamples nC];
    if ~P.save_as_binary
        sig = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/sig'], ...
            dims, maxDims, h5Type, chunkDims, P.deflation);
    else
        % create a binary file where the data is stored
        sig = mysort.ds.binaryFileMatrix(binFile, [1 nC], 'writable', true);
        
        % Save a link to the binary file into /sig:
        binFileName = [name, '.dat'];
        sig_link = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/sig'], [1 length(binFileName) ], [1 length(binFileName) ], 'H5T_C_S1');
        sig_link(1,1:length(binFileName)) = binFileName;
        clear sig_link;
        
        bin_dims = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
        bin_dims(1, :) = [nSamples nC];
        clear bin_dims;
    end
    
    dfun = @(x) int16(x);
    
    chunkSize = 250000;
    nextIdx = 1;
    timer = util.ProcessTimer(ceil(nSamples/chunkSize));
    while lastSample < nSamples
        timer.next();
        timer.showProgress();
        
        thisLastSample = min(nSamples, lastSample+chunkSize);
        
        s = thisLastSample-lastSample;
        X = double(DS(lastSample+1:thisLastSample, :));
        
        %% Check whether there are spikes to be inserted:
        for ui = 1:nU
            st = artificialUnits.spikeTrains{ui, fi};
            sa = artificialUnits.spikeAmplitudes{ui, fi};
            sj = artificialUnits.spikeJitter{ui, fi};
            spikesIdx = find(lastSample < st - cutLeft & st - cutLeft+nTf <= thisLastSample);
            if ~isempty(spikesIdx)
                
                for idx = spikesIdx
                    startFrame = st(idx) - cutLeft - lastSample;
                    lastFrame = st(idx) - cutLeft + nTf - lastSample;
                    X(startFrame+1:lastFrame, :) = X(startFrame+1:lastFrame, :) + sa(idx) * artificialUnits.fp_jittered(:, :, ui, sj(idx));
                end
                
            end
        end
        sig(nextIdx:nextIdx+s-1,:) = dfun(X);
        
        nextIdx = nextIdx+s;
        lastSample = thisLastSample;
    end
    
    % FRAME NUMBERS
    FRAMES = DS.getFrameNumbers();
    
    ffn = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/frame_numbers/first_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
    ffn(1,1) = int32(FRAMES.first_fn);
    clear ffn
    ffn = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/frame_numbers/last_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
    ffn(1,1) = int32(FRAMES.last_fn);
    clear ffn
    ffn = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/frame_numbers/dataDims'], [1 2], [1 2], 'H5T_NATIVE_INT');
    ffn(1,:) = int32(FRAMES.dataDims');
    clear ffn
    ffn = mysort.h5.createVariableAndOrFile(newFileName, [P.sessionName '/frame_numbers/missing_fns'], size(FRAMES.missing_fns'), size(FRAMES.missing_fns'), 'H5T_NATIVE_INT');
    ffn(:,:) = FRAMES.missing_fns';
    clear ffn
    
    disp('Done copying.')
    clear sig
    
    % Set File as being done
    proc(1,1) = int32(0);
    clear proc
end

disp('Generation of files completed')
