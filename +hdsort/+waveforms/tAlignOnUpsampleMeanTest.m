% Templates
T = [];
T(:,1,1) = [0 0  0 -1   -2 -1    4 5 4 3 2 1 0 0];
T(:,1,2) = T(:,1,1)/2;
T(:,1,3) = flipud(T(:,1,1));
T(:,1,4) = [0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0 0 0];
T(:,1,5) = [0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0 0 0 0];
T(:,1,6) = [0 0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0 0];
T(:,1,7) = [0 0 0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0];
T(:,1,7) = [0 0 0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0];
T(:,1,8) = [0 0 0 0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1];
T(:,1,9) = [-1 -1.5 -2 -1.5 .5 1 4 5 1 0 0 0 0 0];

T(:,2,:) = T(:,1,:);

T = T +randn(size(T))*.5;

%%
pd = pdefs();
ff = fullfile(pd.networkTempShare, 'Derk', 'wfs_ForDebuggingAlignment.mat');
wfs = load(ff);
nC = 102;
wfs = hdsort.waveforms.v2t(wfs.wfs, nC);
T = wfs;

[ali tau mMPmask] = hdsort.waveforms.tAlignOnUpsampleMean(T, 'restrictToNMaximalValues', 10, 'initAlignment', []);
% [ali2 tau2 mMPmask] = hdsort.waveforms.tAlignOnUpsampleMean(ali, 'maxIdx', 5);


%%
X = hdsort.waveforms.v4plot(hdsort.waveforms.t2v(T(:,1:2:end,1:4:end)), size(T,2)/2)';
Y = hdsort.waveforms.v4plot(hdsort.waveforms.t2v(ali(:,1:2:end,1:4:end)), size(T,2)/2)';

mysort.plot.figure([1400 800]);
ah = subplot(2,1,1);
plot(X, 'color', [.5 .5 .5]);
hold on
plot(mean(X,2), 'k', 'linewidth', 2);
title('Before Alignment');

ah(2) = subplot(2,1,2);
plot(Y, 'color', [.5 .5 .5]);
hold on
plot(mean(Y,2), 'k', 'linewidth', 2);
title('After Alignment');
linkaxes(ah, 'xy');

% xvsf = hdsort.util.calculateXIvsF(hdsort.waveforms.t2v(T),hdsort.waveforms.t2v(T),1,1);
% mysort.plot.XIvsF(T, T)
% mysort.plot.XIvsF(ali, ali)
