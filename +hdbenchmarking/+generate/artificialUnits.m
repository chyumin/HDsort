function [artificialUnits] = artificialUnits(baseFolder, ...
    sourceFileHandle, swapElectrodePairs, footprintT ...
    parameters)

artificialUnitFile = fullfile(baseFolder, parameters.datasetName, 'artificial_units.mat');

%%
if exist(artificialUnitFile)
    disp('Loading existing artificial units')
    artificialUnits = load(artificialUnitFile);
else
    disp('Create new artificial units')
    
    %% 1. Generate new spiketrains:
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
    
    
    
    
    
    %% Save all the necessary parameters:
    S.parameters.datasetName = parameters.datasetName;
    S.gainMultiplier = double(h5read(DH_original.files.preprocessed{1}, '/Sessions/Session0/filter/gainmultiplier'));
    S.amplitudeCorrectionFactor = parameters.amplitudeCorrectionFactor;
    
    S.swapElectrodePairs = swapElectrodePairs;
    S.electrodeNumbers = ME.electrodeNumbers;
    S.electrodePositions = ME.electrodePositions;
    
    S.parameters.nCellsPerBlock = parameters.nCellsPerBlock;
    S.parameters.spikingRatesHz = parameters.spikingRatesHz;
    S.parameters.sigmaAmplitudes = parameters.sigmaAmplitudes;
    
    S.rawFiles = {};
    for ii = 1:numel(DH_original.files.raw)
        [fpath, name, ext] = fileparts(DH_original.files.raw{ii});
        split_fpath = strsplit(fpath, filesep);
        S.rawFiles{ii} = fullfile('..', split_fpath{end-1}, split_fpath{end}, [name, ext]);
    end
    
    
    
    
    save(artificialUnitFile, '-struct', 'S');
    artificialUnits = S;
    
end