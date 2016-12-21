function [m midx] = vMin(V, nC)
    % computes the channel wise minima and their positions for the
    % vectorized multichannel waveforms in V
    
    % convert to matrix representations, single channel under each other
    M = waveforms.v2m(V,nC);
    % compute channel wise minima and positions
    [m midx] = min(M,[],2);
    % return as vector representation
    m = waveforms.m2v(m, nC);
    midx = waveforms.m2v(midx, nC);