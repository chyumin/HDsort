function normalizedwfs = normalizeEachUnit(wfs)
% Divide each waveform by the amplitude of the negative peak.

negPeakInUnit = min(min(wfs));
n = repmat(negPeakInUnit, size(wfs,1), size(wfs,2), 1);
normalizedwfs = -wfs ./ n;

end