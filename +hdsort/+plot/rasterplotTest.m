% Test Rasterplot:
gdf = [1 20000 1 1;
       1 40000 1 1;
       1 48000 1 1;
       2 24000 1 1;
       2 52000 1 1;
       3 22000 1 1;
       3 34000 1 1;
       3 80000 1 1]

rp = hdsort.plot.Rasterplot(gdf);
footprints = zeros(150, 10, 3);
MultiElectrode = hdsort.file.MultiElectrode([1:10; 1:10]', 1:10)
noiseStd = ones(10, 1);

%%
SPS = hdsort.results.Population('rasterplot_test', ...
    gdf, footprints, MultiElectrode, noiseStd)

rp = hdsort.plot.Rasterplot(SPS)

%%
rp = hdsort.plot.Rasterplot(SPS.Units(2:3))

%%
rp = hdsort.plot.Rasterplot(gdf, ...
    'title', 'Rasterplot',...
    'interval', [10000, 60000], ...
    'unitIds', [1, 3])

