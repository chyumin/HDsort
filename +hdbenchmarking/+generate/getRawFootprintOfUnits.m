function [footprint, P, Q, wfs] = getRawFootprintOfUnits(baseFolder, RAW_list, ST, unitID, varargin)

P.nTf = 150;
P.cutLeft = 50;
P.fir_filterOrder = 110;
P.zeroingThreshold = 4.8;
P.minNChannels = 10;
P.nSubtractMeanChannelIdx = 200;
P = hdsort.util.parseInputs(P, varargin, 'error');

%%
nFiles = numel(RAW_list);
ME = RAW_list(1).MultiElectrode;
nCh = ME.getNElectrodes();
footprint = zeros(P.nTf, nCh);

%% Raw Footprint
bufferFileNewName = fullfile(baseFolder, ['rawFootprintNew' num2str(unitID) '.mat']);
wfsFile = fullfile(baseFolder, ['wfs' num2str(unitID) '.mat']);

try
    load(bufferFileNewName, 'zeroedMeanWFSFiltered', 'meanWFS', ...
        'meanWFSFiltered', 'P', 'Q');
    
    if nargout == 4
        load(wfsFile, 'wfs');
        assert(size(wfs, 3) == Q.nSpikes, 'This must be the same!')
    end
catch
    disp(['Recompute footprint of unit ' num2str(unitID) '...']);
    
    Q.electrodePositions = ME.electrodePositions;
    Q.cutLeftUnfiltered = P.cutLeft + P.fir_filterOrder;
    Q.nTfUnfiltered = P.nTf + 2*P.fir_filterOrder;
    
    if numel(ST) > 5000
        ST = randsample(ST, 5000);
    end
    
    %% Loading cut raw waveforms:
    try
        load(wfsFile);
        
        assert(numel(ST)*0.9 < size(wfs, 3), 'Make sure that at least 90% of all spikes are actually cut out!');
        assert(cutLeft == Q.cutLeftUnfiltered, 'This must be the same!')
        assert(size(wfs,1) == Q.nTfUnfiltered, 'This must be the same!')
    catch
        disp('Cut waveforms...')
        
        wfs_ = []; N = 0; nSpikes = 0;
        for fi = 1:nFiles
            disp(['File ' num2str(fi) ' of ' num2str(nFiles) '...']);
            
            f = RAW_list(fi);
            rel_st = ST(ST >= Q.cutLeftUnfiltered & ST < size(f, 1) - (Q.nTfUnfiltered - Q.cutLeftUnfiltered));
            
            if ~isempty(rel_st)
                wfs_ = [wfs_; f.getWaveform(rel_st, Q.cutLeftUnfiltered, Q.nTfUnfiltered)];
            end
            
            ST = ST - size(f, 1);
            nSpikes = nSpikes + numel(rel_st);
        end
        wfs = mysort.wf.v2t(wfs_, nCh);
        cutLeft = P.cutLeft;
        
        assert(numel(ST)*0.9 < size(wfs, 3), 'Make sure that at least 90% of all spikes are actually cut out!');
        assert(size(wfs, 3) == nSpikes, 'This must be the same!')
        save(wfsFile, 'wfs', 'cutLeft', '-v7.3')
        
    end
    Q.nSpikes = size(wfs, 3);
    
    
    meanWFS = mean(wfs, 3);
    
    %% Bandpass filter the mean waveforms:
    [meanWFSFiltered, filterProperties] = hdbenchmarking.generate.filterFootprints(meanWFS, 'fir_filterOrder', P.fir_filterOrder);
    Q.filterProperties = filterProperties;
    
    %% Subtract the mean of the waveforms on each channel and assertain a smooth onset and offset:
    [zeroedMeanWFSFiltered, zeroingProperties] = hdbenchmarking.generate.assertainZeroOffset(meanWFSFiltered);
    Q.zeroingProperties = zeroingProperties;
    
    %%
    save(bufferFileNewName, 'zeroedMeanWFSFiltered',  'P', 'Q', ...
        'meanWFSFiltered', 'meanWFS', '-v7.3');
    
end

%% Return the footprint:
footprint = zeroedMeanWFSFiltered;

end