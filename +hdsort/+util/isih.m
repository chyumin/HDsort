function [isih times_ms P] = isih(hdsort.spiketrain., varargin)
if ~iscell(hdsort.spiketrain.)
%    if size(hdsort.spiketrain., 1) == 1
%        hdsort.spiketrain. = hdsort.spiketrain.(:);
%    end
%    assert( size(hdsort.spiketrain., 2) == 1, 'Input must be a single column vector!')
    if size(hdsort.spiketrain., 2) > 1
        hdsort.spiketrain. = hdsort.spiketrain.fromGdf(hdsort.spiketrain.);
    else
        [isih times_ms P] = hdsort.util.isih({hdsort.spiketrain.}, varargin{:});
        return
    end
end

P.binSize_ms = 1.0;
P.maxlag_ms = 20.0; % ms
P.Fs = 20000;
P = hdsort.util.parseInputs(P, varargin, 'error');

assert(iscell(hdsort.spiketrain.), 'Input variable must be a cell array!')

maxlag = P.maxlag_ms / 1000.0 * P.Fs;
binSize = P.binSize_ms / 1000.0 * P.Fs;
nBins = round( maxlag / binSize);
bins = linspace(0, maxlag, nBins);

times_ms = bins / P.Fs * 1000.0;

nTrials = length(hdsort.spiketrain.);
for ii = 1:nTrials
    isis = diff(hdsort.spiketrain.{ii});
    isis = isis(isis <= maxlag);
    isih(ii, :) = hist(isis, bins);
end

end