function [V rangeidx] = v4hdsort.plot.V, nC)
    V = hdsort.waveforms.m2v([hdsort.waveforms.v2m(V, nC) nan(size(V,1)*nC,1)], nC);
    if nargout == 2
        rangeidx = hdsort.waveforms.vSubChannelIdx(V, nC);
    end