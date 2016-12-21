function [mmi, mma, mi_ch, ma_ch] = tGlobalMinMax(T)
    % Returns the global minima and maxima for each template in T with the
    % channel numbers
    % Input: 
    %    T  - tensor (time x channels x templates)
    % Output:
    %    mmi - minima for each template
    %    mma - maxima
    %    mi_ch - channel index for that minimum
    %    ma_ch - same for max
    
    [mi, ma, mi_idx, ma_idx] = waveforms.tMinMaxPerTemplate(T);
    [mma ma_ch] = max(ma');
    [mmi mi_ch] = min(mi');
    
    
    
    