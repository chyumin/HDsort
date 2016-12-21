
function [tau A] = alignWaveformsUpsampleMin(S, nC, varargin)
    warning('This function is depricated! Use mysort.wf.* instead!');
    [tau A] = mysort.hdsort.util.alignWaveformsUpsampleMax(-S, nC, varargin{:});
    A = -A;
