function [maxChannels, maxVal] = maxChannel(fp)
if size(fp, 2) == 1
    fp_ = [fp fp*0];
    [maxChannels, maxVal] = waveforms.maxChannel(fp_);
    maxChannels = 0*maxChannels + 1;
    return
end

[m m_] = max(abs(fp));
[m2 maxChannels] = max(m);
maxChannels = maxChannels(:);
maxVal = squeeze(m2);