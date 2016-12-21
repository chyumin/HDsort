T = [];
T(1:10,1,1) = [1:10];
T(1:10,2,1) = [1:10]+1;

T(1:10,1,2) = -[1:10];
T(1:10,2,2) = -[1:10]+1;

T(1:10,1,3) = -[1:10]+11;
T(1:10,2,3) = -[1:10]+11+1;

gdf = [1 10
       2 30
       3 50
       1 100
       2 105
       3 110];
cutLeft = 3;
[x y] = hdsort.waveforms.templateSpikeTrainData(T, gdf, cutLeft); 

figure;
hdsort.plot.x,y(1,:), 'b')
hold on
hdsort.plot.x,y(2,:), 'g')