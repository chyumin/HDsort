function [newFileList, artificialUnitFile] = artificialUnits(...
    datasetName, ...
    baseFolder, ...
    gdf, ... % matrix containing the unitID in the first row and the spike-times in the second
    amplitudes, ... % Give the amplitudes in noiseStd that the estimated footprints should have
    estimated_footprints, ... % the characteristic waveforms of each unit on all channels
    RAW_list, .... % handle on raw file (to extract the footprints)
    PRE, ... % handle on pre-filtered file
    swapElectrodePairs, ... % from the function that creates the new MultiElectrode
    blockIdx, ...
    varargin) % assignment to which block an electrode in MultiElectrode belongs to

%% Parameters:
parameters.randomSeed = 1;
parameters.nCellsPerBlock = 10;
parameters.spikingRatesHz = 5;
parameters.sigmaAmplitudes = 0.1;
parameters.refractoryPeriod_s = 0.0015;
parameters.jitterFactor = 10;

parameters.unitselection.footprintSelectionCriterium = 'random'; % | 'targetamplitude' --> 1.5 mean amplitude of all units
parameters.unitselection.forbiddenUnits = []; % Specify a list of units that should be excluded from the possible candidates

%parameters.footprints.amplitudeCorrectionFactor = 0.5; % footprintSelectionCriterium == 'random' --> 0.5 ; 1.0 otherwise
parameters.footprints.nTf = 150;
parameters.footprints.cutLeft = 50;
parameters.footprints.zeroingThreshold = 4.8;
parameters.footprints.minNChannels = 10;
parameters.footprints.nSubtractMeanChannelIdx = 200;
parameters.footprints.filter.hpf = 300;
parameters.footprints.filter.lpf = 6000;
parameters.footprints.filter.fir_filterOrder = 110;
parameters.footprints.filter.frameRate = 20000;
parameters.footprints.filter.doZeroPad = false;

parameters.datafile.save_as_binary = true;
parameters.datafile.deleteOldFilesIfNecessary = false;
parameters = hdsort.util.parseInputs(parameters, varargin, 'error');

%% 0. Prepare generation:
unitIDs = unique(gdf(:,1));
original_spiketrains = hdsort.spiketrain.fromGdf(gdf);
MultiElectrode = PRE.MultiElectrode;

mkdir(fullfile(baseFolder, datasetName));
rng(parameters.randomSeed);


%% 1. Select units that you want to use as templates for the artificial units:
selectedUnitsFile = fullfile(baseFolder, datasetName, 'selected_units.mat');
if exist(selectedUnitsFile)
    disp('Loading existing selected units...')
    selectedUnits = load(selectedUnitsFile);
else
    selectedUnits = hdbenchmarking.generate.selectUnits( gdf, estimated_footprints, ...
        parameters.nCellsPerBlock, blockIdx, parameters.unitselection)
    save(selectedUnitsFile, '-struct', 'selectedUnits')
end

%% 2. Swap footprints:
footprintsFile = fullfile(baseFolder, datasetName, 'artificial_footprints.mat');
if exist(footprintsFile)
    disp('Loading existing artificial footprints...')
    FP = load(footprintsFile);
    disp('Done.')
else
    %%
    disp('Create new artificial footprints...')
    bufferFolder = fullfile(baseFolder, 'buffer');
    mkdir(bufferFolder)
    
    FP.swapElectrodePairs = swapElectrodePairs;
    
    FP.originalCellIdxBlock0 = selectedUnits.originalCellIdxBlock0;
    FP.originalCellIdxBlock1 = selectedUnits.originalCellIdxBlock1;
    FP.unitIdx = [FP.originalCellIdxBlock0 FP.originalCellIdxBlock1];
    
    %% Compute full footprints for the selected units:
    mf = RAW_list(1).getMissingFrameNumbers();
    [fp_, P_, Q_] = hdbenchmarking.generate.getRawFootprintOfUnits(bufferFolder, RAW_list, [mf.first + 1000], 1, parameters.footprints);
    footprints_unswapped = zeros(size(fp_, 1), size(fp_, 2), 2*parameters.nCellsPerBlock);
    AMPs = amplitudes([FP.originalCellIdxBlock0; FP.originalCellIdxBlock1]);
    
    
    indices = [FP.originalCellIdxBlock0; FP.originalCellIdxBlock1]';
    parfor ii = 1:length(indices) 
        unitIdx = indices(ii);
        hdbenchmarking.generate.getRawFootprintOfUnits(bufferFolder, RAW_list, ...
            original_spiketrains{unitIdx}, unitIDs(unitIdx), parameters.footprints);
    end
    for ii = 1:parameters.nCellsPerBlock
        unitIdx = FP.originalCellIdxBlock0(ii);
        [fp_, footprintP(ii), footprintQ(ii)] = ...
            hdbenchmarking.generate.getRawFootprintOfUnits(bufferFolder, RAW_list, ...
            original_spiketrains{unitIdx}, unitIDs(unitIdx), parameters.footprints);
        
        %% Adjust the amplitudes:
        footprints_unswapped(:, :, ii) = fp_; % * parameters.footprints.amplitudeCorrectionFactor;
    end
    for jj = parameters.nCellsPerBlock+1:2*parameters.nCellsPerBlock
        
        unitIdx = FP.originalCellIdxBlock1(jj-parameters.nCellsPerBlock);
        [fp_, footprintP(jj), footprintQ(jj)] = ...
            hdbenchmarking.generate.getRawFootprintOfUnits(bufferFolder, RAW_list, ...
            original_spiketrains{unitIdx}, unitIDs(unitIdx), parameters.footprints);
        
        %% Adjust the amplitudes:
        footprints_unswapped(:, :, jj) = fp_; %* parameters.footprints.amplitudeCorrectionFactor;
    end
    
    %% Normalize:
    footprints_unswapped_normalized_ = hdsort.waveforms.normalizeEachUnit(footprints_unswapped);
    AMP_rep = repmat( reshape(AMPs, 1, 1, length(AMPs)), ...
        size(footprints_unswapped_normalized_, 1), ...
        size(footprints_unswapped_normalized_, 2), 1); 
    footprints_unswapped_normalized = AMP_rep .* footprints_unswapped_normalized_;
    FP.noiseStd = mean(PRE.noiseStd);
    
    footprints_unswapped = FP.noiseStd * footprints_unswapped_normalized;
    FP.footprintP = footprintP;
    FP.footprintQ = footprintQ;
    FP.cutLeft = footprintP(1).cutLeft;
    
    %% Swap the footprints:
    swappedFP = hdbenchmarking.generate.swapFootprint(footprints_unswapped, MultiElectrode, FP.swapElectrodePairs);
    [nTf, nCh, nU] = size(swappedFP);
    FP.footprints = swappedFP;
    
    %% Create jittered footprints:
    fpUP = hdsort.waveforms.tResample(FP.footprints, parameters.jitterFactor, 1, 1);
    FP.fp_jittered = zeros(nTf, nCh, nU, parameters.jitterFactor);
    for j = 1:parameters.jitterFactor
        FP.fp_jittered(:, :, :, j) = fpUP(j:parameters.jitterFactor:end, :,:);
    end
    
    save(footprintsFile, '-struct', 'FP')
end

%% 3. Generate new spiketrains:
spikeTrainsFile = fullfile(baseFolder, datasetName, 'artificial_spiketrains.mat');
if exist(spikeTrainsFile)
    disp('Loading existing spike trains...')
    ST = load(spikeTrainsFile);
    disp('Done.')
else
    disp('Create new spike trains...')
    
    nArtificialUnits = 2 * parameters.nCellsPerBlock;
    if numel(parameters.spikingRatesHz) == 1
        parameters.spikingRatesHz = repmat(parameters.spikingRatesHz, 1, nArtificialUnits);
    else
        parameters.spikingRatesHz = parameters.spikingRatesHz( randperm(numel( parameters.spikingRatesHz)) );
    end
    assert( numel(parameters.spikingRatesHz) == nArtificialUnits, 'Give a spiking rate for each artificial unit!')
    
    fn = double(PRE.getFrameNumbers());
    samplingRate = PRE.getSampleRate();
    
    ST = hdbenchmarking.generate.spikeTrains(samplingRate, fn(1), fn(end), ...
        parameters.spikingRatesHz, parameters.sigmaAmplitudes, ...
        parameters.refractoryPeriod_s, parameters.jitterFactor);
    
    assert( ~any(~ST.spikeCount(:)), 'There must be at least one spike per file!');
    
    save(spikeTrainsFile, '-struct', 'ST')
    disp('Done.')
end

%% 4. Create the artificial units and save them in one file:
artificialUnitFile = fullfile(baseFolder, datasetName, 'artificial_units.mat');
if exist(artificialUnitFile)
    disp('Loading existing artificial units')
    artificialUnits = load(artificialUnitFile);
    disp('Done.')
else
    disp('Create new artificial units')

    %% Save all the necessary parameters:
    artificialUnits.datasetName = datasetName;
    artificialUnits.parameters = parameters;
    artificialUnits.FP = FP;
    artificialUnits.ST = ST;
    artificialUnits.MultiElectrode = MultiElectrode; 
    
    artificialUnits.rawFiles = {};
    for ii = 1:numel(RAW_list)
        [fpath, name, ext] = fileparts(RAW_list(ii).fileName);
        split_fpath = strsplit(fpath, filesep);
        artificialUnits.rawFiles{ii} = fullfile('..', split_fpath{end-1}, split_fpath{end}, [name, ext]);
    end
    
    save(artificialUnitFile, '-struct', 'artificialUnits');
    disp('Done.')
end

%% 5. Create the new files by inserting the artificial units into the pre-processed files:
outputLocation = fullfile(baseFolder, datasetName);
newFileList = hdbenchmarking.generate.dataFiles(PRE, outputLocation, ...
                                    artificialUnitFile, parameters.datafile)

end                             