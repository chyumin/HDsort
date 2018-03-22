function ST = spikeTrains(samplingRate, firstFrames, lastFrames,...
                    spikingRates_Hz, sigmaAmplitudes, refractoryPeriod_s, jitterFactor)
% Input:
% spikingRateHz - For Poisson random spiketime generator
% firstFrames / lastFrames - Framenumbers for each file for absolute spiketimes
% nCells - Number of cells
% sigmaAmplitudes - Standard deviation of spiking amplitudes
% fileName - specify where the variable are saved to
% samplingRate - sampling rate
% refractoryPeriod_s - refractory period in s
% 
% Output:
% spikeTrains - nCell x nFiles cell array containing vectors of spiketimes
% spikeTrainsAbsTime - same as 'spikeTrains' but times are relative to files
% spikeAmplitudes - cell array specifying a relative amplitude for each spike

nCells = numel(spikingRates_Hz);

% sigmaAmplitudes: variance of the normalized amplitude (around 1)
assert(length(firstFrames) == length(lastFrames), 'Error: lastFrames must be same length as firstFrames!')

durations = [];
spikeTrains = {};
spikeTrainsAbsTime = {};
spikeAmplitudes = {};
spikeJitter = {};
for jj = 1:nCells
    for ii = 1:length(lastFrames)
        durations(ii) = lastFrames(ii)-firstFrames(ii)+1;
        spikeTrains{jj, ii} = round( poisson_spike_train(spikingRates_Hz(jj)/samplingRate, durations(ii)) );
        
        if ~isempty(refractoryPeriod_s)
            idx = diff(spikeTrains{jj, ii}) < refractoryPeriod_s*samplingRate;
            spikeTrains{jj, ii}(idx) = [];
        end
        
        spikeTrainsAbsTime{jj, ii} = spikeTrains{jj, ii} + firstFrames(ii); % in frames relative to files
        spikeAmplitudes{jj, ii} = randn(length(spikeTrains{jj, ii}),1)* sigmaAmplitudes + 1;
        spikeJitter{jj, ii} = ceil(rand(length(spikeTrains{jj, ii}), 1)*jitterFactor);
    end
end
info = 'First Dimension: Cells - Second Dimension: Files';

ST = struct;
ST.spikeTrains = spikeTrains;
ST.spikeTrainsAbsTime = spikeTrainsAbsTime;
ST.spikeAmplitudes = spikeAmplitudes;
ST.spikeJitter = spikeJitter;

ST.spikeCount = cellfun(@numel, spikeTrains);

ST.sigmaAmplitudes = sigmaAmplitudes;
ST.spikingRates_Hz = spikingRates_Hz;
ST.jitterFactor = jitterFactor;

ST.durations = durations;

ST.firstFrames = firstFrames;
ST.lastFrames = lastFrames;

ST.nCells = nCells;
ST.samplingRate = samplingRate;
ST.refractoryPeriod_s = refractoryPeriod_s;

ST.info = info;

end