%%
X = rand(1, 100);
Y = rand(1, 100);
groupIdx = randi(5, [1 100]);

%%
gs = myplot.Gscatter(X,Y,groupIdx);

%%
gs = myplot.Gscatter(X,Y,groupIdx, 'markerType', 'o');
%%
gs = myplot.Gscatter(X,Y,groupIdx, 'color', [0 0 0]);

