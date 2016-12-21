function [A, tau] = vAlignOnUpsampleMean(S, nC, varargin)
    A = waveforms.v2t(S,nC);
    [A, tau] = waveforms.tAlignOnUpsampleMean(A, varargin{:});
    A = waveforms.t2v(A);