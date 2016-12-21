
%%
a = repmat(1:10, 1, 5);
[idx, y, ny, P] = hdsort.util.samplePopulation(a, 3)

%% With a cell array as input:
c = {};
c{1} = 1:10;
c{2} = 30:50;
c{3} = (20:40)'+(60:80)'*i;
[idx, y, ny, P] = hdsort.util.samplePopulation(c, 3)

%% Throw error when input type not consistent:
c{4} = (1:10)==3
[idx, y, ny, P] = hdsort.util.samplePopulation(c, 3)


%%
d = logical(zeros(1,100));
d(2) = true;
d(4) = true;
[idx, y, ny, P] = hdsort.util.samplePopulation(d, 3)

idx(logical(y))