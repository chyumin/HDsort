baseFolder = './Tutorial';

%% Inputs:
% Sorting results:
gdf = R.gdf_merged;
original_footprint =  R.T_merged;
unitIDs = unique(gdf(:,1));
original_spiketrains = hdsort.spiketrain.fromGdf(gdf);

rawFileNames
PRE     % preprocessed file handle
MultiElectrode = PRE.MultiElectrode;
swapElectrodePairs % from the function that creates the new MultiElectrode
blockIdx % assignment to which block an electrode in MultiElectrode belongs to

%% Parameters:
parameters.datasetName = 'amplitude_sweep01';
parameters.randomSeed = 1;
parameters.amplitudeCorrectionFactor = 0.5; % selectSameAmplitudeFootprints == false --> 0.5 ; 1.0 otherwise
parameters.nCellsPerBlock = 10;
parameters.spikingRatesHz = 5;
parameters.sigmaAmplitudes = 0.1;
parameters.refractoryPeriod_s = 0.0015;
parameters.jitterFactor = 10;

parameters.unitselection.footprintSelectionCriterium = 'random'; % | 'targetamplitude' --> 1.5 mean amplitude of all units
parameters.unitselection.forbiddenUnits = [1, 3, 5]; % Specify a list of units that should be excluded from the possible candidates

%% 0. Prepare generation:
mkdir(fullfile(baseFolder, parameters.datasetName));
rng(parameters.randomSeed)


%% 1. Select units that you want to use as templates for the artificial units:
selectedUnitsFile = fullfile(baseFolder, parameters.datasetName, 'selected_units.mat');
if exist(selectedUnitsFile)
    disp('Loading existing selected units...')
    selectedUnits = load(selectedUnitsFile);
else
    selectedUnits = hdbenchmarking.generate.selectUnits( gdf, original_footprint, ...
        parameters.nCellsPerBlock, blockIdx, parameters.unitselection)
    save(selectedUnitsFile, '-struct', 'selectedUnits')
end

%% 2. Swap footprints:
footprintsFile = fullfile(baseFolder, parameters.datasetName, 'artificial_footprints.mat');
if exist(footprintsFile)
    disp('Loading existing artificial footprints...')
    FP = load(footprintsFile);
else
    
    disp('Create new artificial footprints...')
    FP.gainMultiplier = PRE.getGainMultiplier();
    FP.swapElectrodePairs = swapElectrodePairs;
    
    FP.originalCellIdxBlock0 = selectedUnits.originalCellIdxBlock0;
    FP.originalCellIdxBlock1 = selectedUnits.originalCellIdxBlock1;
    FP.unitIdx = [FP.originalCellIdxBlock0 FP.originalCellIdxBlock1];
    
    %% Compute full footprints for the selected units:
    [fp_, P_, Q_] = hdbenchmarking.generate.getRawFootprintOfUnits(baseFolder, RAW, [1000], 1);
    footprints_unswapped = zeros(size(fp_, 1), size(fp_, 2), 2*parameters.nCellsPerBlock);
    
    for ii = 1:parameters.nCellsPerBlock
        unitIdx = FP.originalCellIdxBlock0(ii);
        [fp_, footprintP(ii), footprintQ(ii)] = ...
            hdbenchmarking.generate.getRawFootprintOfUnits(baseFolder, RAW, original_spiketrains{unitIdx}, unitIDs(unitIdx));
        
        %% Adjust the amplitudes:
        footprints_unswapped(:, :, ii) = fp_ * parameters.amplitudeCorrectionFactor;
    end
    for jj = parameters.nCellsPerBlock+1:2*parameters.nCellsPerBlock
        
        unitIdx = FP.originalCellIdxBlock0(jj);
        [fp_, footprintP(jj), footprintQ(jj)] = ...
            hdbenchmarking.generate.getRawFootprintOfUnits(baseFolder, RAW, original_spiketrains{unitIdx}, unitIDs(unitIdx));
        
        %% Adjust the amplitudes:
        footprints_unswapped(:, :, jj) = fp_ * parameters.amplitudeCorrectionFactor;
    end
    footprints_unswapped = FP.gainMultiplier * footprints_unswapped;
    FP.footprintP = footprintP;
    FP.footprintQ = footprintQ;
    FP.cutLeft = footprintP(1).cutLeft;
    
    %% Swap the footprints:
    [swappedFP, blockIdx] = hdbenchmarking.generate.swapFootprint(footprints_unswapped, MultiElectrode, FP.swapElectrodePairs);
    [nTf, nCh, nU] = size(swappedFP);
    FP.footprints = swappedFP;
    
    %% Create jittered footprints:
    fpUP = mysort.wf.tResample(FP.footprints, FP.jitterFactor, 1, 1);
    FP.fp_jittered = zeros(nTf, nCh, nU, FP.jitterFactor);
    for j = 1:FP.jitterFactor
        FP.fp_jittered(:, :, :, j) = fpUP(j:FP.jitterFactor:end, :,:);
    end
    
    save(footprintsFile, '-struct', 'FP')
end

%% 3. Generate new spiketrains:
spikeTrainsFile = fullfile(baseFolder, parameters.datasetName, 'artificial_spiketrains.mat');
if exist(spikeTrainsFile)
    disp('Loading existing spike trains...')
    ST = load(spikeTrainsFile);
else
    disp('Create new spike trains...')
    
    nArtificialUnits = 2 * parameters.nCellsPerBlock;
    if numel(parameters.spikingRatesHz) == 1
        parameters.spikingRatesHz = repmat(parameters.spikingRatesHz, 1, nArtificialUnits);
    else
        parameters.spikingRatesHz = parameters.spikingRatesHz( randperm(numel( parameters.spikingRatesHz)) );
    end
    assert( numel(parameters.spikingRatesHz) == nArtificialUnits, 'Give a spiking rate for each artificial unit!')
    
    fn = sourceFileHandle.getFrameNumbers();
    samplingRate = sourceFileHandle.getSampleRate();
    
    ST = hdbenchmarking.generate.spiketrains(samplingRate, fn.first_fn, fn.last_fn, ...
        parameters.spikingRatesHz, parameters.sigmaAmplitudes, ...
        parameters.refractoryPeriod_s, parameters.jitterFactor);
    
    assert( ~any(~ST.spikeCount(:)), 'There must be at least one spike per file!');
    
    save(spikeTrainsFile, '-struct', 'ST')
end

%% 4. Create the artificial units and save them in one file:

artificialUnitFile = fullfile(baseFolder, parameters.datasetName, 'artificial_units.mat');
if exist(artificialUnitFile)
    disp('Loading existing artificial units')
    artificialUnits = load(artificialUnitFile);
else
    disp('Create new artificial units')

    %% Save all the necessary parameters:
    artificialUnits.parameters = parameters;
    artificialUnits.FP = FP;
    artificialUnits.ST = ST;
    artificialUnits.MultiElectrode = MultiElectrode; 
    
    artificialUnits.rawFiles = {};
    for ii = 1:numel(rawFileNames)
        [fpath, name, ext] = fileparts(rawFileNames);
        split_fpath = strsplit(fpath, filesep);
        artificialUnits.rawFiles{ii} = fullfile('..', split_fpath{end-1}, split_fpath{end}, [name, ext]);
    end
    
    save(artificialUnitFile, '-struct', 'artificialUnits');
end

%%