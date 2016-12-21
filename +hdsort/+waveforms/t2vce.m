function vce = t2vce(t)
    [Tf, nC, nT] = size(t);
    
    if isempty(t)
        vce = [];
        return
    end
    
    vce = hdsort.waveforms.vte2vce(hdsort.waveforms.t2v(t), nC);