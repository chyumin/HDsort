% Defines
nC = 1000;
L  = 10000;

% create simulated data
X = randn(nC,L);

% make variance of channel 2 bigger
X(2,:) = 2*X(2,:);
% correlate channel 3 and 4
X(4,:) = X(4,:) + X(3,:);

% build the covest:
invertFullMatrix = 1;
fc_Full = mysort.hdsort.noise.FakeCovest(X,invertFullMatrix);
invertFullMatrix = 0;
fc_Diag = mysort.hdsort.noise.FakeCovest(X,invertFullMatrix);

% prewhiten a vector:
Tf = 40;
x = [ones(1, Tf*nC); randn(1, Tf*nC)];
y_full = fc_Full.invMul(x);
y_diag = fc_Diag.invMul(x);


%%
figure;
subhdsort.plot.2,3,1)
imagesc(fc_Full.xcovs(1:10,1:10))
title('Top Left Part of C')
ylabel('Full Covariance');
subhdsort.plot.2,3,2)
imagesc(fc_Full.iC(1:10,1:10))
title('Top Left Part of inv(C)')
subhdsort.plot.2,3,3)
hdsort.plot.y_full(:,1:200)');
title('Prewhitened Vector of ones')
subhdsort.plot.2,3,1+3)
imagesc(diag(diag(fc_Diag.xcovs(1:10,1:10))))
ylabel('Only Diag Covariance');
subhdsort.plot.2,3,2+3)
imagesc(fc_Diag.iC(1:10,1:10))
subhdsort.plot.2,3,3+3)
hdsort.plot.y_diag(:,1:200)');