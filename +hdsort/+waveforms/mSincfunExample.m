%%
x = [0:.1:8*pi];
y = [2*sin(x);
     2*cos(x)];
sf = hdsort.waveforms.mSincfun(x, y);
xi1 = 0:.05:2*pi;
yi1 = sf(xi1);

xi2 = 0:.2:2*pi;
yi2 = sf(xi2);

xi3 = 8+(0:.2:4*pi);
yi3 = sf(xi3);

figure;
subhdsort.plot.2,1,1)
hdsort.plot.x, y(1,:), '.-b');
hold on
hdsort.plot.xi1, yi1(1,:), '.-r');
hdsort.plot.xi2, yi2(1,:), '.-c');
hdsort.plot.xi3, yi3(1,:), '.-g');
subhdsort.plot.2,1,2)
hdsort.plot.x, y(2,:), '.-b');
hold on
hdsort.plot.xi1, yi1(2,:), '.-r');
hdsort.plot.xi2, yi2(2,:), '.-c');
hdsort.plot.xi3, yi3(2,:), '.-g');


%%
x = [0:.1:8*pi];
y = [x;
     -x];
sf = hdsort.waveforms.mSincfun(x, y);
xi1 = 0:.05:2*pi;
yi1 = sf(xi1);

xi2 = 0:.2:2*pi;
yi2 = sf(xi2);

xi3 = 8+(0:.2:4*pi);
yi3 = sf(xi3);

figure;
subhdsort.plot.2,1,1)
hdsort.plot.x, y(1,:), '.-b');
hold on
hdsort.plot.xi1, yi1(1,:), '.-r');
hdsort.plot.xi2, yi2(1,:), '.-c');
hdsort.plot.xi3, yi3(1,:), '.-g');
subhdsort.plot.2,1,2)
hdsort.plot.x, y(2,:), '.-b');
hold on
hdsort.plot.xi1, yi1(2,:), '.-r');
hdsort.plot.xi2, yi2(2,:), '.-c');
hdsort.plot.xi3, yi3(2,:), '.-g');

%%
figure
x = -10:.05:10;
y = mysort.hdsort.util.sinc0(x);
hdsort.plot.x,y)
