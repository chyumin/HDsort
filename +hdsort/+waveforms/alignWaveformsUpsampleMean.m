
function [tau A] = alignWaveformsUpsampleMean(S, nC, varargin)
    warning('This function is depricated! Use mysort.wf.* instead!');
    %P.maxIdx = [];
    P.UP = 3;
    P.nIter = 3;
    P.downsample = 1;
    P = mysort.hdsort.util.parseInputs(P, varargin, 'error');

    Sup = mysort.hdsort.util.resampleTensor(mysort.hdsort.util.m2t(S, nC), P.UP,1);
    Sup = mysort.hdsort.util.t2m(Sup);
    [tau Sup] = mysort.hdsort.util.alignWaveforms(Sup, nC, 'nIter', P.nIter);%, 'maxIdx', P.UP*P.maxIdx);
    if P.downsample
        tau = tau/P.UP;
        A = mysort.hdsort.util.resampleTensor(mysort.hdsort.util.m2t(Sup, nC), 1, P.UP);
        A = mysort.hdsort.util.t2m(A);
    else
        A = Sup;
    end