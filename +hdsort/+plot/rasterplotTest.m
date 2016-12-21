% Test Rasterplot:
gdf = [1 20000;
       1 40000;
       1 48000;
       2 24000;
       2 52000;
       3 22000;
       3 34000;
       3 80000]

rp = myplot.Rasterplot(gdf);

%%
SPS = lsa.SpikeSorting(gdf)
rp = myplot.Rasterplot(SPS.Units)

%%
rp = myplot.Rasterplot(gdf, ...
    'title', 'Rasterplot',...
    'interval', [10000, 60000], ...
    'unitIds', [1, 3])

