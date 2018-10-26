function [binnedST, edges, edges_seconds] = binSpikeTrains(ST, window, sampling_rate, interval)
% Return for each 'spiketrain' (vector of spikeedges_seconds) the spikecount of each
% bin defined by 'interval'
% The intervals are [interval(i), interval(i+1)[

if ~iscell(ST)
    ST = {ST};
end

assert(size(ST, 2) == 1, 'Input spiketrains must be in a Nx1 cell array!')

if nargin < 4
    allst = cat(1, ST{:});
    allst = allst(:);
    interval = [min(allst) max(allst)];
end

dt = window*sampling_rate;
edges = interval(1):dt:interval(2);
binnedST = hdsort.spiketrain.toBins(ST, edges);
edges_seconds = edges / sampling_rate;

end