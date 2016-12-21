function [V rangeidx] = v4plot(V, nC)
    V = waveforms.m2v([waveforms.v2m(V, nC) nan(size(V,1)*nC,1)], nC);
    if nargout == 2
        rangeidx = waveforms.vSubChannelIdx(V, nC);
    end