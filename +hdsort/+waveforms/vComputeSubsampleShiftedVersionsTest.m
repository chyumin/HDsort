
T = [0 1 2 3 0 0 1 2 3 0
     0 3 2 1 0 0 3 2 1 0];
nC = 2;

    
vT = hdsort.waveforms.vComputeSubsampleShiftedVersions(T, nC, 3);

figure;
subhdsort.plot.2,1,1);
hdsort.plot.vT(:,:,1)');
subhdsort.plot.2,1,2);
hold on
for i=1:size(vT,3)
    hdsort.plot.vT(:,:,i)', 'color', mysort.hdsort.plot.vectorColor(i))
end

%%
xrange = 0:.4:2*pi;
T = [sin(xrange) sin(xrange)
     cos(xrange*2) cos(xrange*2)];
nC = 2;
    
vT = hdsort.waveforms.vComputeSubsampleShiftedVersions(T, nC, 20);

figure;
subhdsort.plot.2,1,1);
hdsort.plot.vT(:,:,1)');
subhdsort.plot.2,1,2);
hold on
for i=1:size(vT,3)
    hdsort.plot.vT(:,:,i)', 'color', mysort.hdsort.plot.vectorColor(i))
end
