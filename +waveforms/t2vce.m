function vce = t2vce(t)
    [Tf, nC, nT] = size(t);
    
    if isempty(t)
        vce = [];
        return
    end
    
    vce = waveforms.vte2vce(waveforms.t2v(t), nC);