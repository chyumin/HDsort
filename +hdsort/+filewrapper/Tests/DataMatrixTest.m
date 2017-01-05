%% DataMatrixTest

M = randn(10000, 100);
MM = hdsort.filewrapper.DataMatrix(M);

%%
MM.restrictToChannels();
wfs0 = MM.getWaveform(100, 20, 50, 1:3);
st0 = MM.detectSpikes();
x0 = MM.getData(1:10,1:3);

MM.restrictToChannels(1:2);
wfs1 = MM.getWaveform(100, 20, 50);
st1 = MM.detectSpikes();
x1 = MM.getData(1:10,:);

MM.restrictToChannels();
MM.restrictToChannels(2:3);
wfs2 = MM.getWaveform(100, 20, 50);
st2 = MM.detectSpikes();
x2 = MM.getData(1:10,:);

%%
assert(~any( wfs1(51:100) - wfs2(1:50) ) , 'error!')
assert(~any( wfs0(51:100) - wfs1(51:100) ) , 'error!')

%%
assert(~any( st1{2} - st2{1}) , 'error!')
assert(~any( st1{1} - st0{1}) , 'error!')

%%
assert(~any( x1(:, 2) -  x2(:, 1) ) , 'error!')
assert(~any( x1(:, 1) -  x0(:, 1) ) , 'error!')




