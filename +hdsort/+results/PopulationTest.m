% Test Rasterplot:
gdf = [1 20000 1 2;
       1 40000 1 1;
       1 48000 1 1;
       2 24000 1 10;
       2 52000 1 1;
       3 22000 1 4;
       3 34000 1 8;
       3 80000 1 1;
       3 80100 1 3;
       3 80200 1 6;
       3 80300 1 4;
       3 80400 1 2];

footprints = zeros(150, 10, 3);
MultiElectrode = hdsort.file.MultiElectrode([1:10; 1:10]', 1:10)
noiseStd = (1:10) + 0.1;

SP = hdsort.results.Population('rasterplot_test', ...
    gdf, footprints, MultiElectrode, noiseStd)


gdf1 = SP.getGdf();
gdf1(:,1:2)

%% Merge units 1 and 2 into unit 4
U4 = SP.mergeUnits([1,2]);
gdf2 = SP.getGdf();
gdf2(:,1:2)

%% Split unit 3 into unit 5 and 6

Usplit = SP.splitUnit(3, [1,1,0,2,0,1,1]);

gdf3 = SP.getGdf();
gdf3(:,1:2)

%%
SP.splitmerge;

%%
U7 = SP.mergeUnits([1, 2]);

gdf4 = SP.getGdf();
gdf4(:,1:2)

%%
SP.revertMerge(U7)
gdf5 = SP.getGdf();
gdf5(:,1:2)

assert( isequal(gdf5, gdf3), '!')

%%

SP.revertSplit(Usplit)
gdf6 = SP.getGdf();
gdf6(:,1:2)

assert( isequal(gdf2, gdf6), '!')

%%
SP.revertMerge(U4)
gdf7 = SP.getGdf();
gdf7(:,1:2)

assert( isequal(gdf7, gdf1), '!')

% ---- Back to beginning ---

%% 
Usplit2 = SP.splitUnit(3, [1,1,0,0,1,1,1]);

gdf8 = SP.getGdf();
gdf8(:,1:2)

Usplit3 = SP.splitUnit(5, [1,1,0,0,1]);
gdf9 = SP.getGdf();
gdf9(:,1:2)

%%
SP.revertSplit(Usplit3)
gdf10 = SP.getGdf();
gdf10(:,1:2)


