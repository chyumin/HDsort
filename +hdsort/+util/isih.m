function [isih times_ms P] = isih(hdsort_spiketrain, varargin)
if ~iscell(hdsort_spiketrain)
%    if size(hdsort_spiketrain, 1) == 1
%        hdsort_spiketrain = hdsort_spiketrain(:);
%    end
%    assert( size(hdsort_spiketrain, 2) == 1, 'Input must be a single column vector!')
    if size(hdsort_spiketrain, 2) > 1
        hdsort_spiketrain = hdsort.spiketrain.fromGdf(hdsort_spiketrain);
    else
        [isih times_ms P] = hdsort.util.isih({hdsort_spiketrain}, varargin{:});
        return
    end
end

P.binSize_ms = 1.0;
P.maxlag_ms = 20.0; % ms
P.Fs = 20000;
P = hdsort.util.parseInputs(P, varargin, 'error');

assert(iscell(hdsort_spiketrain), 'Input variable must be a cell array!')

maxlag = P.maxlag_ms / 1000.0 * P.Fs;
binSize = P.binSize_ms / 1000.0 * P.Fs;
nBins = round( maxlag / binSize);
bins = linspace(0, maxlag, nBins);

times_ms = bins / P.Fs * 1000.0;

nTrials = length(hdsort_spiketrain);
for ii = 1:nTrials
    isis = diff(hdsort_spiketrain{ii});
    isis = isis(isis <= maxlag);
    isih(ii, :) = hist(isis, bins);
end

end