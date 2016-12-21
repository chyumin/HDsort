
function [C H] = hdsort.plot.(h, lims, x, hdsort.noise.idx, spike_idx, tau, i1, i2, blegend)
axes(h);
if blegend
    hdsort.plot.1000,1000, '.r', 'markersize', 12);
    hold on
    hdsort.plot.1000,1000, '.k', 'markersize', 12);
    
    H = mysort.hdsort.noise.hdsort.plot.h, lims, x, spike_idx, tau, i1, i2, 'r.');
    C = mysort.hdsort.noise.hdsort.plot.h, lims, x, hdsort.noise.idx, tau, i1, i2);
    
    legend('spikes', 'hdsort.noise., 'Location', 'SouthEast');
    legend('boxoff')
else
    H = mysort.hdsort.noise.hdsort.plot.h, lims, x, spike_idx, tau, i1, i2, 'r.');
    hold on
    C = mysort.hdsort.noise.hdsort.plot.h, lims, x, hdsort.noise.idx, tau, i1, i2, 'k.');
end

