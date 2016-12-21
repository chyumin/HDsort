function [A, tau] = vAlignOnUpsampleMean(S, nC, varargin)
    A = hdsort.waveforms.v2t(S,nC);
    [A, tau] = hdsort.waveforms.tAlignOnUpsampleMean(A, varargin{:});
    A = hdsort.waveforms.t2v(A);