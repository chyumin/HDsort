

%% With no holes in the config:
[X1, Y1] = meshgrid(1:10, 1:10);
[X2, Y2] = meshgrid(31:40, 31:40);
X = [X1(:), Y1(:); X2(:), Y2(:)] * 17.5;
X = X(randperm(size(X,1)),:);
ME = hdsort.file.MultiElectrode(X, 1:size(X,1));

[newME, swapPairsEl, blockIdx] = hdbenchmarking.generate.overlappingBlocks(ME);

newME.plotConfig();
figure; gscatter(newME.electrodePositions(:,1), newME.electrodePositions(:,2), blockIdx) 

%% With holes:
[X1, Y1] = meshgrid(1:10, 1:10);
[X2, Y2] = meshgrid(31:40, 31:40);
X = [X1(:), Y1(:); X2(:), Y2(:)] * 17.5;

remove_els = [1, 11, 22, 33, 102, 112, 122, 108, 117, 127, 137, 147];
X(remove_els,:) = [];
X = X(randperm(size(X,1)),:);
ME = hdsort.file.MultiElectrode(X, 1:size(X,1));

ME.plotConfig()

[newME, swapPairsEl, blockIdx] = hdbenchmarking.generate.overlappingBlocks(ME)

newME.plotConfig();
figure; gscatter(newME.electrodePositions(:,1), newME.electrodePositions(:,2), blockIdx) 


