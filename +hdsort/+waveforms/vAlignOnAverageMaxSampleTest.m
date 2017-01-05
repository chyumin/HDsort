    nS = 100;
    nC = 2;

    X = repmat([5:-.5:0 1:10 11 10:-1:1 0:.5:5], nS, nC)-5;
    maxIdx = 22;
    X = X + 1*randn(size(X));
    shifts = [-4:2:4];
    stau = repmat(shifts', nS/length(shifts),1);
    figure; plot(X(1,:));
    X = hdsort.util.shiftMCRows(X,stau,nC,1);
    hold on; plot(X(1,:), 'r');
%     Y1 = hdsort.util.shiftMCRows(X,-stau,nC,1);
%     Y2 = hdsort.util.shiftRowsInterpolated(X, -stau , nC);
    debug = 1;
%     mysortx.plot.spikes(X,'nC',nC);
%     title('Not aligned Spikes'); 
%     mysortx.plot.spikes(Y1,'nC',nC); title('Back perf. Aligned Spikes'); 
%     mysortx.plot.spikes(Y2,'nC',nC); title('Back perf. interp Aligned Spikes'); 
%     [tau Y] = hdsort.util.alignWaveformsOnMax(X, 2, 'debug',1);
% 
%     fprintf('Shift Error after Iteration %d: %d\n', i, sum(abs(-tau-stau)));
%     mysortx.plot.spikes(Y,'nC',nC);
%     title('Aligned Spikes');
    
    [tau Y] = hdsort.waveforms.vAlignOnAverageMaxSample(X, nC, 'maxIdx', 22, 'nMaxChannelsForWeighting', 5);
    fprintf('Shift Error: %.4f\n', sum(abs(-tau-stau)));
    figure; plot(stau, -tau, '.');
    xlabel('real tau');
    ylabel('estimated tau');
    mysortx.plot.spikes(Y,'nC',nC);
    title('Aligned Spikes');    
%     [tau Y] = hdsort.util.alignWaveformsOnMax(X, 2, 'debug',1, 'maxIdx', 30);
%     fprintf('Shift Error after Iteration %d: %d\n', i, sum(abs(-tau-stau)));
%     mysortx.plot.spikes(Y,'nC',nC);
%     title('Aligned Spikes');        

