function V = vSubsel(V, nC, idx)
    % subselects the time index into multi channel hdsort.waveforms.on every
    % channel if the hdsort.waveforms.are stored in concatenated form as rows of
    % matrix V
    if isempty(V)
        return
    end
    vidx = hdsort.waveforms.vSubIdx(V,nC,idx);
    V = V(:,vidx);