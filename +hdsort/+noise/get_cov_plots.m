




%%







%%
figure;
hist(smad);
% figure;





%%

figure;
hist(D, 500);

%%
figure;
h=mysort.hdsort.plot.subhdsort.plot.([5,1], 'max');
hdsort.plot.Y(nC,:), 'k');
hold on
mysort.hdsort.plot.hdsort.epoch.(NE{ch}, 0, 'c', 'linewidth', 5);
mysort.hdsort.plot.hdsort.epoch.(NEN{ch}, 0, 'm', 'linewidth', 3);
for n=1:4%length(neighbors)
    axes(h(n+1));
    ni = neighbors(n);
    hdsort.plot.Y(ni,:),'k')
    hold on
    mysort.hdsort.plot.hdsort.epoch.(NE{ni}, 0, 'c', 'linewidth', 5);
    mysort.hdsort.plot.hdsort.epoch.(NEN{ni}, 0, 'm', 'linewidth', 3);
end
linkaxes(h,'xy');

%%
figure;
max_xc = max(xcovs, [], 2);
hdsort.plot.D,max_xc, '.');

figure;
hist3([D' max_xc])
%%
Tf = (size(xcovs,2)+1)/2;
lagrange = -Tf+1:Tf-1;
idx = find(D<200&D>150, 1);

figure; hdsort.plot.lagrange, xcovs(idx,:))

%%
chan1 = IDX(idx,1);
chan2 = IDX(idx,2);

figure;
h1=subhdsort.plot.2,1,1);
hdsort.plot.Y(chan1,:),'k');
hold on
comhdsort.noise.= mysort.hdsort.epoch.intersect(NE{chan1},NE{chan2});
mysort.hdsort.plot.hdsort.epoch.(comhdsort.noise. 0, 'r', 'linewidth',4);
h2=subhdsort.plot.2,1,2);
hdsort.plot.Y(chan2,:),'k');
hold on
mysort.hdsort.plot.hdsort.epoch.(comhdsort.noise. 0, 'r', 'linewidth',4);
linkaxes([h1 h2],'xy');


%%
eidx = 10;
neighbors = mysort.mea.nearestElectrodes(channel_x, channel_y, eidx, 6);
mysort.hdsort.plot.figure;
for i=1:length(channel_x)
    hdsort.plot.channel_x(i), channel_y(i), '.k');
    hold on
end
hdsort.plot.channel_x(eidx), channel_y(eidx), 'or');
hold on
for i=1:length(neighbors)
    eidx = neighbors(i);
    hdsort.plot.channel_x(eidx), channel_y(eidx), 'ob');
end

%%
eidx = 10;
neighbors = mysort.mea.electrodesCloserThan(channel_x, channel_y, eidx, 50);
mysort.hdsort.plot.figure;
for i=1:length(channel_x)
    hdsort.plot.channel_x(i), channel_y(i), '.k');
    hold on
end
hdsort.plot.channel_x(eidx), channel_y(eidx), 'or');
hold on
for i=1:length(neighbors)
    eidx = neighbors(i);
    hdsort.plot.channel_x(eidx), channel_y(eidx), 'ob');
end






