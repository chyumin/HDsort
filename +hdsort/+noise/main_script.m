

% %% Buffer short
% X = X(:, 1:400000);
% save('C:\LocalData\Michele\marching_square_buffer_short.mat', 'X','srate','channel_x', 'channel_y');


%% Alternative ALLE
% H = mysort.hdsort.util.hdf5recursiveLoad('C:\LocalData\Alle\AnalyseAlleCompleteHDF5\Preprocessing\buf0003.h5');
% X1 = H.RESi.r(1,1).X;
% X2 = H.RESi.r(1,2).X;
% save('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1', 'X2');

%% Load data
% % 
% clear all
% close all
% load('C:\LocalData\Michele\marching_square_buffer_short.mat', 'X','srate','channel_x', 'channel_y');
% k1 = 1;
% k2 = 4;

% Alternative
% load('C:\LocalData\Alle\temp_prep_buf0001h5.mat', 'X');
% k1 = 1;
% k2 = 2;

% Alternative
% clear all
% close all
% load('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1');
% X = X1; clear X1;
% k1 = 3;
% k2 = 4;

% Alternative
% clear all
% close all
% load('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1', 'X2');
% shift = 1000;
% X = [X1(3,1:end-shift); X2(3,1+shift:end)]; clear X1 X2;
% k1 = 1;
% k2 = 2;

%% Get sigma
Tf = 100;
thr1 = 4.5;
[smad1, spikehdsort.epoch.1] = ana.douglas.estimateSigma(X(k1,:), Tf, thr1);
hdsort.noise.dsort.epoch1 = mysort.hdsort.epoch.flip(spikehdsort.epoch.1, size(X,2));
[smad2, spikehdsort.epoch.2] = ana.douglas.estimateSigma(X(k2,:), Tf, thr1);
hdsort.noise.dsort.epoch2 = mysort.hdsort.epoch.flip(spikehdsort.epoch.2, size(X,2));

spikehdsort.epoch.1_idx = mysort.hdsort.epoch.toIdx(spikehdsort.epoch.1);
spikehdsort.epoch.2_idx = mysort.hdsort.epoch.toIdx(spikehdsort.epoch.2);
hdsort.noise.dsort.epoch1_idx = mysort.hdsort.epoch.toIdx(hdsort.noise.dsort.epoch1);
hdsort.noise.dsort.epoch2_idx = mysort.hdsort.epoch.toIdx(hdsort.noise.dsort.epoch2);
commonhdsort.noise.dsort.epoch = mysort.hdsort.epoch.intersect(hdsort.noise.dsort.epoch1,hdsort.noise.dsort.epoch2);
commonhdsort.noise.dsort.epoch_idx = mysort.hdsort.epoch.toIdx(commonhdsort.noise.dsort.epoch);
unionspikehdsort.epoch.  = mysort.hdsort.epoch.merge( [spikehdsort.epoch.1; spikehdsort.epoch.2]);
unionspikehdsort.epoch._idx = mysort.hdsort.epoch.toIdx(unionspikehdsort.epoch.);
%% show test of data
fig = mysort.hdsort.plot.figure();
h1 = subhdsort.plot.2,1,1);
hdsort.plot.X(k1,:), 'k'); hold on

thr2 = thr1;
hdsort.plot.[0 size(X,2)], thr2*[smad1 smad1], ':g', 'linewidth', 2)
hdsort.plot.[0 size(X,2)],-thr2*[smad1 smad1], ':g', 'linewidth', 2)
mysort.hdsort.plot.hdsort.epoch.(unionspikehdsort.epoch., 0, 'm', 'linewidth', 6);
mysort.hdsort.plot.hdsort.epoch.(spikehdsort.epoch.1, 0, 'c', 'linewidth', 4);
mysort.hdsort.plot.hdsort.epoch.(hdsort.noise.dsort.epoch1, 0, 'r', 'linewidth', 4);
mysort.hdsort.plot.hdsort.epoch.(commonhdsort.noise.dsort.epoch, 0, 'k', 'linewidth', 2);


h2 = subhdsort.plot.2,1,2);
hdsort.plot.X(k2,:), 'k'); hold on

hdsort.plot.[0 size(X,2)], thr2*[smad2 smad2], ':g', 'linewidth', 2)
hdsort.plot.[0 size(X,2)],-thr2*[smad2 smad2], ':g', 'linewidth', 2)
mysort.hdsort.plot.hdsort.epoch.(unionspikehdsort.epoch., 0, 'm', 'linewidth', 6);
mysort.hdsort.plot.hdsort.epoch.(spikehdsort.epoch.2, 0, 'c', 'linewidth', 4);
mysort.hdsort.plot.hdsort.epoch.(hdsort.noise.dsort.epoch2, 0, 'r', 'linewidth', 4);
mysort.hdsort.plot.hdsort.epoch.(commonhdsort.noise.dsort.epoch, 0, 'g', 'linewidth', 2);

linkaxes([h1 h2],'xy');
%% Estimate autocov
maxLag = 50;
L = size(X,2);
r1 = xcorr(X(k1,:), maxLag, 'biased');
r2 = xcorr(X(k2,:), maxLag, 'biased');
xd = xcorr(X(k1,:), X(k2,:), maxLag, 'biased');

[c1 Ln1] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X(k1,:), hdsort.noise.dsort.epoch1, maxLag, maxLag);
[c2 Ln2] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X(k2,:), hdsort.noise.dsort.epoch2, maxLag, maxLag);
[xc Lxn] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X([k1 k2],:),commonhdsort.noise.dsort.epoch, maxLag, maxLag);

[sp1 Lsp1] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X(k1,:),spikehdsort.epoch.1, maxLag, maxLag);
[sp2 Lsp2] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X(k2,:),spikehdsort.epoch.2, maxLag, maxLag);
[xsp Lxsp] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(X([k1 k2],:),unionspikehdsort.epoch., maxLag, maxLag);


%% Plot covs
fig = mysort.hdsort.plot.figure();
xr = -maxLag:maxLag;
hdsort.plot.xr, r1, 'k:', 'linewidth', 2); hold on
hdsort.plot.xr, r2, 'c:', 'linewidth', 2);
hdsort.plot.xr, xc, 'm:', 'linewidth', 2);

hdsort.plot.xr, c1, 'k', 'linewidth', 2); 
hdsort.plot.xr, c2, 'c', 'linewidth', 2);
hdsort.plot.xr, xd, 'm', 'linewidth', 2);

hdsort.plot.xr, sp1, 'k--', 'linewidth', 2); 
hdsort.plot.xr, sp2, 'c--', 'linewidth', 2);
hdsort.plot.xr, xsp, 'm--', 'linewidth', 2);
legend('datacov1', 'datacov2', 'data_xcov', 'hdsort.noise.ov1', 'hdsort.noise.ov2', 'hdsort.noise.xcov', 'spcov1', 'spcov2', 'sp_xcov');

%% Test relationships
fig = mysort.hdsort.plot.figure();
hdsort.plot.xr, c1*Ln1, 'b');
hold on
hdsort.plot.xr, sp1*Lsp1, 'r');
hdsort.plot.xr, r1*L, 'g');
hdsort.plot.xr, c1*Ln1 + sp1*Lsp1, 'k');
residual = c1*Ln1 + sp1*Lsp1 - r1*L;
hdsort.plot.xr, residual, 'k:');

legend('N', 'SP', 'R', 'N+SP', 'err');


%%
tau = 1;
x = [X(k1,:)/smad1; X(k2,:)/smad2];
lims = 15*[-1 1];
%% Plot instantanen template space between hdsort.noise.in both channels
fig = mysort.hdsort.plot.figure('width', 740, 'height', 800);
h = mysort.hdsort.plot.subhdsort.plot.([2,2]);

[Ck1k1_1 Hk1k1_1] = mysort.hdsort.noise.hdsort.plot.(h(1), lims, x, hdsort.noise.dsort.epoch1_idx, spikehdsort.epoch.1_idx, tau, 1, 1, true);
[Ck2k2_1 Hk2k2_1] = mysort.hdsort.noise.hdsort.plot.(h(3), lims, x, hdsort.noise.dsort.epoch2_idx, spikehdsort.epoch.2_idx, tau, 2, 2, false);
[Ck1k2_0 Hk1k2_0] = mysort.hdsort.noise.hdsort.plot.(h(2), lims, x, commonhdsort.noise.dsort.epoch_idx, unionspikehdsort.epoch._idx,   0, 1, 2, false);
[Ck1k2_1 Hk1k2_1] = mysort.hdsort.noise.hdsort.plot.(h(4), lims, x, commonhdsort.noise.dsort.epoch_idx, unionspikehdsort.epoch._idx, tau, 1, 2, false);
 


%% OLD PLOT
% fig = mysort.hdsort.plot.figure('width', 1000, 'height', 600);
% h1 = subhdsort.plot.2,4,1);
% Ck1k1_1 = mysort.hdsort.noise.hdsort.plot.h1, lims, x, hdsort.noise.dsort.epoch1_idx, tau, 1, 1, 'k.', 'hdsort.noise.);
% h2 = subhdsort.plot.2,4,2);
% Hk1k1_1 = mysort.hdsort.noise.hdsort.plot.h2, lims, x, spikehdsort.epoch.1_idx, tau, 1, 1, 'r.', 'spikes');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subhdsort.plot.2,4,5);
% Ck2k2_1 = mysort.hdsort.noise.hdsort.plot.h1, lims, x, hdsort.noise.dsort.epoch2_idx, tau, 2, 2);
% h2 = subhdsort.plot.2,4,6);
% Hk2k2_1 = mysort.hdsort.noise.hdsort.plot.h2, lims, x, spikehdsort.epoch.2_idx, tau, 2, 2, 'r.');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subhdsort.plot.2,4,3);
% Ck1k2_0 = mysort.hdsort.noise.hdsort.plot.h1, lims, x, commonhdsort.noise.dsort.epoch_idx, 0, 1, 2);
% h2 = subhdsort.plot.2,4,4);
% Hk1k2_0 = mysort.hdsort.noise.hdsort.plot.h2, lims, x, unionspikehdsort.epoch._idx, 0, 1, 2, 'r.');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subhdsort.plot.2,4,7);
% Ck1k2_1 = mysort.hdsort.noise.hdsort.plot.h1, lims, x, commonhdsort.noise.dsort.epoch_idx, tau, 1, 2);
% h2 = subhdsort.plot.2,4,8);
% Hk1k2_1 = mysort.hdsort.noise.hdsort.plot.h2, lims, x, unionspikehdsort.epoch._idx, tau, 1, 2, 'r.');
% linkaxes([h1 h2], 'xy');

