
%% ------------------------------------------------------------------------
maxRate = 50;
rate_time = 1:0.05:30;
rate = maxRate * 0.5 * (sin(5*rate_time)+1);
%%
train = {};
for ii = 1:100
    train{ii} = hdsort.util.poissonSpikeTrain(rate, rate_time);
end

%%
train = [];
for ii = 1:100
    train_ = hdsort.util.poissonSpikeTrain(rate, rate_time);
    train = [train train_];
end
train = sort(train);
figure; hist(train, 200)
%%

myhdsort.plot.Rasterplot(train)

%%
figure; hist(diff(train{1}), 10000)

%% ------------------------------------------------------------------------
L = 100000; %ms
freq = 1./linspace(1, 100, L);
freq(freq< 0.75) = 0;
fftx = freq .* (1 + i*2*pi*randn(size(freq))) ;
intensity = real(ifft(fftx));
intensity = (intensity - min(intensity)) / (max(intensity)- min(intensity));

rate_time = 0.001*(1:length(intensity));
rate = intensity*50;

%%
train = {};
for ii = 1:1000
    train{ii} = hdsort.util.poissonSpikeTrain(rate, rate_time);
end

%%
train_ = [];
for ii = 1:numel(train)
    train_ = [train_, train{ii}];
end

figure;
a(1) = subplot.2,1,1)
hist(train_, rate_time)
a(2) = subplot.2,1,2)
hdsort.plot.rate_time, rate*1000)
linkaxes(a, 'x')

