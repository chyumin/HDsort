function [templates, nSpikesPerTemplate, P] = templateEstimation(DS, gdf, varargin)

P.Tf = 150;
P.cutLeft = 50;
P.maxSpikes = 100;
P.wfsFile = '';
P.nSpikesChunk = 1000;
P.bufferFolder = '';
P.roundSpikeTimes = true;
P = hdsort.util.parseInputs(P, varargin, 'error');

if isempty(gdf)
    templates = [];
    nSpikesPerTemplate = [];
    return;
end

if P.roundSpikeTimes
    % This avoids an annoying warning
    gdf(:,2) = round(gdf(:,2));
end

assert(size(DS,1) >= max(gdf(:,2)), 'Spiketimes outside of range!')

if ~isempty(P.wfsFile)
    if exist(P.wfsFile, 'file') ~= 2
        gdf_reduced = reduceGDF(gdf, P.maxSpikes);
        
        %% Prepare the waveform file:
        nC = size(DS, 2);
        samplesPerSecond = DS.samplesPerSecond;
        cutLeft = P.cutLeft;
        Tf = P.Tf;
        W = hdsort.waveforms.WaveFormFile(P.wfsFile, 'cutLeft', cutLeft, ...
            'Tf', Tf, 'nC', nC, 'samplesPerSecond', samplesPerSecond);
        
        %% Cut the waveforms:
        disp('Cut waveforms...')
        chunker = hdsort.util.Chunker(size(gdf_reduced, 1), 'chunkSize', P.nSpikesChunk, 'progressDisplay', 'console');
        while chunker.hasNextChunk()
            chunk = chunker.getNextChunk();
            chunked_gdf = gdf_reduced(chunk(1):chunk(2), :);
            wfs_ = DS.getWaveform(chunked_gdf(:,2), cutLeft, Tf);
            W.addWaveforms(wfs_, chunked_gdf);
        end
    else
        %% Load precut waveforms:
        W = hdsort.waveforms.WaveFormFile(P.wfsFile);
        gdf_reduced = W.getGdf();
        
        Tf = W.getTf();
        nC = W.getNChannels();
        cutLeft = W.getCutLeft();
        
        assert( cutLeft == P.cutLeft, 'Waveforms not saved in the correct format!')
        assert( Tf == P.Tf, 'Waveforms not saved in the correct format!')
    end
    
    % #################### Actually estimate templates ####################
    checkFormat(gdf_reduced, gdf, P.maxSpikes)
    tic
    [templates, nSpikesPerTemplate] = estimateTemplates(W, gdf_reduced, Tf, nC);
    toc
else
    [gdf_reduced, unitIDs] = reduceGDF(gdf, P.maxSpikes);
    checkFormat(gdf_reduced, gdf, P.maxSpikes);
    
    nC = size(DS, 2);
    Tf = P.Tf;
    cutLeft = P.cutLeft;
    maxBytes = 12e9;
    
    if size(gdf_reduced, 1) < maxBytes / (8*Tf*nC)
        %% Cut the waveforms in a simple way when there are not too many:
        disp('Cut waveforms (simple)...')
        if ~isempty(gdf_reduced)
            fprintf('Cutting Spikes on all channels to estimate full templates...\n');
            t = tic;
            W = DS.getWaveform(gdf_reduced(:,2), cutLeft, Tf);
            t = toc(t);
            fprintf('Cutting Spikes done. (%.2fs)\n', t);
        else
            W = [];
        end
        
        % #################### Actually estimate templates ####################
        tic
        [templates, nSpikesPerTemplate] = estimateTemplates(W, gdf_reduced, Tf, nC);
        toc
    else
        %% Cut the waveforms and average them to templates chunk-wise if there are too many:
        assert(~isempty(P.bufferFolder) && exist(P.bufferFolder, 'dir')==7, 'If you want to estimate so many templates at the same time, you must provide a buffer folder!')
        bufferFile = fullfile(P.bufferFolder, 'templates_buffer.mat');
        
        templates = []; nSpikesPerTemplate = []; nSavedChunks = 0;
        
        nUnitsMaxPerChunk = floor((maxBytes /(8*Tf*nC)) / P.maxSpikes);
        N = 0;
        chunker_over_units = hdsort.util.Chunker(numel(unitIDs), 'chunkSize', nUnitsMaxPerChunk, 'progressDisplay', 'console');
        while chunker_over_units.hasNextChunk()
            unitsIdx = chunker_over_units.getNextChunk();
            N = N + 1;
            
            %% Skip cutting if the buffer file has already got content
            if nSavedChunks >= N
                continue;
            end
            if exist(bufferFile, 'file')==2
                load(bufferFile);
            end
            if nSavedChunks >= N
                continue;
            end
            
            unitIDs_chunk = unitIDs(unitsIdx(1):unitsIdx(2));
            gdf_reduced_over_units = gdf_reduced( ismember(gdf_reduced(:,1),unitIDs_chunk) , :);
            
            if ~isempty(gdf_reduced_over_units)
                fprintf('Cutting Spikes on all channels to estimate full templates...\n');
                t = tic;
                W = DS.getWaveform(gdf_reduced_over_units(:,2), cutLeft, Tf);
                t = toc(t);
                fprintf('Cutting Spikes done. (%.2fs)\n', t);
            else
                W = [];
            end
            
            % #################### Actually estimate templates ####################
            tic
            [templates_, nSpikesPerTemplate_] = estimateTemplates(W, gdf_reduced_over_units, Tf, nC);
            toc
            assert( size(templates_, 3) == numel(unitIDs_chunk), 'There must be the same number of templates as units!')
            
            % Save them to the output variables:
            templates = cat(3, templates, templates_);
            nSpikesPerTemplate = cat(1, nSpikesPerTemplate, nSpikesPerTemplate_);
            nSavedChunks = nSavedChunks + 1;
            save(bufferFile, 'templates', 'nSpikesPerTemplate', 'nSavedChunks');
            clear W templates_ nSpikesPerTemplate_
        end
        load(bufferFile);
        assert( size(templates, 3) == numel(unitIDs), 'There must be the same number of templates as units!')
    end
end

% ---------------------------------------------------------------------
    function checkFormat(gdf_reduced, gdf, maxSpikes)
        [uIDs_, nuIDs_] = hdsort.util.unique(gdf_reduced(:,1));
        assert( ~any( nuIDs_ > maxSpikes ), 'There are too many waveforms for some units!')
        assert( isequal(uIDs_(:), unique(gdf(:,1))), 'Some units seem to be missing!')
    end

% ---------------------------------------------------------------------
    function [templates, nSpikesPerTemplate] = estimateTemplates(W, gdf_reduced, Tf, nC)
        %% Estimate the templates:
        uIDs = unique(gdf_reduced(:,1));
        nUnits = numel(uIDs);
        
        templates = zeros(Tf, nC, nUnits);
        nSpikesPerTemplate = zeros(nUnits, 1);
        for ii = 1:nUnits
            if nUnits > 20
                disp([num2str(ii) ' of ' num2str(nUnits) ' computed...'])
            end
            if isempty(gdf_reduced) continue; end
            idx = find(gdf_reduced(:,1) == uIDs(ii));
            templates(:,:,ii) = hdsort.waveforms.v2t( median(W(idx, :),1), nC);
            nSpikesPerTemplate(ii) = numel(idx);
        end
        disp('Template estimation finished.')
    end
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
    function [gdf_reduced, uIDs] = reduceGDF(gdf, maxSpikes)
        %% Prepare the spiketimes:
        uIDs = unique(gdf(:,1));
        nUnits = numel(uIDs);
        
        gdf_reduced = [];
        for ii = 1:nUnits
            idx = find(gdf(:,1) == uIDs(ii));
            nSpikes = numel(idx);
            if maxSpikes < nSpikes
                idx = randsample(idx, maxSpikes);
            end
            gdf_reduced = [gdf_reduced; gdf(idx, :)];
        end
        [~, sidx] = sort(gdf_reduced(:,2));
        gdf_reduced = gdf_reduced(sidx, :);
    end
% ---------------------------------------------------------------------

end