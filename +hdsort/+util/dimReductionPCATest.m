% template1 = sin(0:.01:2*pi);
% template1 = template1/max(abs((template1)));
D = load('dimReductionPCATestTemplate.mat');
templates = D.templates;
templates = resample(waveforms.vSubChanSel(templates, D.nC, 4)', 11, 1)'; 
template1 = templates(1,:)/max(abs((templates(1,:))));
template2 = templates(2,:)/max(abs((templates(2,:))));
nC = 1;
Tf_org = length(template1);
down = round(Tf_org/55);
T = [template1' template2'];
T1      = repmat([template1 zeros(1,down)]', 1, down);
T1([1:down end-down+1:end],:) = [];
Tshift1 = toeplitz([template1 zeros(1,down)], [template1(1) zeros(1,down-1)]);
Tshift1([1:down end-down+1:end],:) = [];
T2      = repmat([template2 zeros(1,down)]', 1, down);
T2([1:down end-down+1:end],:) = [];
Tshift2 = toeplitz([template2 zeros(1,down)], [template2(1) zeros(1,down-1)]);
Tshift2([1:down end-down+1:end],:) = [];
TT = {[T1 T2], [Tshift1 Tshift2]};

for it = 1:length(TT)
    myT = TT{it};
    Tdown = resample(myT, 1, down);
    Tdownali      = waveforms.vAlignOnMax(Tdown', nC)';
    % make more versions
    Tdown = repmat(Tdown, 1, 100);
    N = .1*randn(size(Tdown'));
    [TdownaliNoise taus] = waveforms.vAlignOnMax(Tdown'+ N, nC);
    TdownaliNoise = TdownaliNoise';
    TdownaliNoise([1:max(1,max(taus)) end+min(0,min(taus))+1:end],:) = [];
    
%     [taus TdownaliIntNoise] = waveforms.vAlignOnAverageMaxSample(Tdown'+ N, 1);
    [TdownaliIntNoise taus] = waveforms.vAlignOnUpsampleMean(Tdown'+ N, nC, 'maxIter', 3);
    TdownaliIntNoise = TdownaliIntNoise';
    TdownaliIntNoise([1:max(1,ceil(max(taus))) end+min(0,min(floor(taus)))+1:end],:) = [];

    %%
    figure;
    nR = 5; nCol = 3;
    subplot.nR, nCol,1+nCol*0)
    hdsort.plot.T); axis tight;
    title('Template only')
    subplot.nR, nCol,1+nCol*1)
    hdsort.plot.Tdown); axis tight;
    title('Template With Shifts Downsampled')
    subplot.nR, nCol,1+nCol*2)
    hdsort.plot.Tdownali); axis tight;
    title('Template With Shifts Downsampled and Aligned')
    subplot.nR, nCol,1+nCol*3)
    hdsort.plot.TdownaliNoise); axis tight;
    title('Template With Shifts Downsampled, added Noise and then Aligned')
    subplot.nR, nCol,1+nCol*4)
    hdsort.plot.TdownaliIntNoise); axis tight;
    title('Template With Shifts Downsampled, added Noise and then Aligned and Int')
    %%
    dim = 2;
    [XFet1 pcs1 T1] = hdsort.util.dimReductionPCA(T', dim, [], [], 1);
    [XFet1down pcs1down T1down] = hdsort.util.dimReductionPCA(Tdown', dim, [], [], 1);
    [XFet1downali pcs1downali T1downali] = hdsort.util.dimReductionPCA(Tdownali', dim, [], [], 1);
    [XFet1downaliNoise pcs1downaliNoise T1downaliNoise] = hdsort.util.dimReductionPCA(TdownaliNoise', dim, [], [], 1);
    [XFet1downaliIntNoise pcs1downaliIntNoise T1downaliIntNoise] = hdsort.util.dimReductionPCA(TdownaliIntNoise', dim, [], [], 1);

    subplot.nR, nCol,2+nCol*0)
    hdsort.plot.XFet1(:,1), XFet1(:,2), '.');
    subplot.nR, nCol,2+nCol*1)
    hdsort.plot.XFet1down(:,1), XFet1down(:,2), '.');
    subplot.nR, nCol,2+nCol*2)
    hdsort.plot.XFet1downali(:,1), XFet1downali(:,2), '.');
    subplot.nR, nCol,2+nCol*3)
    hdsort.plot.XFet1downaliNoise(:,1), XFet1downaliNoise(:,2), '.');
    subplot.nR, nCol,2+nCol*4)
    hdsort.plot.XFet1downaliIntNoise(:,1), XFet1downaliIntNoise(:,2), '.');
    
    upsample_factor = 11;
    [XFet3 pcs3 T3] = hdsort.util.dimReductionPCA(T', dim, [], [], upsample_factor);
    [XFet3down pcs3down T3down] = hdsort.util.dimReductionPCA(Tdown', dim, [], [], upsample_factor);
    [XFet3downali pcs3downali T3downali] = hdsort.util.dimReductionPCA(Tdownali', dim, [], [], upsample_factor);
    [XFet3downaliNoise pcs3downaliNoise T3downaliNoise T3shiftPCs] = hdsort.util.dimReductionPCA(TdownaliNoise', dim, [], [], upsample_factor);
    [XFet3downaliIntNoise pcs3downaliIntNoise T3downaliIntNoise] = hdsort.util.dimReductionPCA(TdownaliIntNoise', dim, [], [], upsample_factor);

    subplot.nR, nCol,3+nCol*0)
    hdsort.plot.XFet3(:,1), XFet3(:,2), '.');
    subplot.nR, nCol,3+nCol*1)
    hdsort.plot.XFet3down(:,1), XFet3down(:,2), '.');
    subplot.nR, nCol,3+nCol*2)
    hdsort.plot.XFet3downali(:,1), XFet3downali(:,2), '.');
    subplot.nR, nCol,3+nCol*3)
    hdsort.plot.XFet3downaliNoise(:,1), XFet3downaliNoise(:,2), '.');
    subplot.nR, nCol,3+nCol*4)
    hdsort.plot.XFet3downaliIntNoise(:,1), XFet3downaliIntNoise(:,2), '.');    
end
%%
figure
subplot.3,1,1)
hdsort.plot.pcs1downaliIntNoise(:,1:2));
title('Correct PCs');
subplot.3,1,2)
hdsort.plot.pcs3downaliNoise(:,1:2));
title('Approx PCs');
subplot.3,1,3)
hdsort.plot.squeeze(T3shiftPCs(:,1:2,1)));
title('Shifted Approx PCs');