%%
sp1 = hdsort.plot.Subplots([2 3], 'title', 'subplot1')
sp2 = hdsort.plot.Subplots([2 3], 'title', 'subplot2')
%
gs1 = hdsort.plot.Gscatter(rand(1, 100),rand(1, 100),randi(5, [1 100]), 'title', 'gs1');
gs2 = hdsort.plot.Gscatter(rand(1, 100),rand(1, 100),randi(5, [1 100]), 'title', 'gs2');

%
sp1.add(gs1, 1);
sp1.add(gs2, 2);

%
sp1.add(hdsort.plot.Gscatter(rand(1, 100),rand(1, 100),randi(5, [1 100]), 'title', 'gs3'))
%
sp2.add(gs1);
sp1.show


%% ---
sp1 = hdsort.plot.Subplots([2 3], 'title', 'subplot1')
sp1.add(hdsort.plot.Gscatter(rand(1, 100),rand(1, 100),randi(5, [1 100]), 'title', 'gs1'), 1);
sp1.add(hdsort.plot.Gscatter(rand(1, 100),rand(1, 100),randi(5, [1 100]), 'title', 'gs2'), 2);
