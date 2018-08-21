
%%
rate = 100; % Hz
duration = 60*100; % s
train = hdbenchmarking.generate.poissonSpikeTrain(rate, duration);

%
if 1
    refractoryPeriod_s = 0.003;
    idx = [false, diff(train) < refractoryPeriod_s];
    while sum(idx)
        train(idx) = [];
        idx = [false, diff(train) < refractoryPeriod_s];
    end
end

%%
[isih, times_ms, P] = hdsort.util.isih({train*20000}, 'maxlag_ms', 100, 'binSize_ms', 0.5);
figure; bar(times_ms, isih/sum(isih), 'histc')
xlabel('isih [ms]')

%%
x = 0:10; % ms
figure; plot(x, 1./factorial(x), '-o')


%%

%%

% interspike_time probability:
rate = 50;
probability = 0.0:0.01:1.0; % s
interspike_time = -log(probability)/rate
%interspike_time20 = -log(probability)/20
%interspike_time50 = -log(probability)/50

probability( interspike_time < refractoryPeriod_s) = 0;

figure; hold on;
plot(interspike_time, probability)
%plot(interspike_time50, probability)
%plot(interspike_time20, probability)
%plot(interspike_time10*10, probability)
xlabel('interspike time difference [ms]')
ylabel('probability')
legend([ num2str(rate) 'Hz'])

%plot([refractoryPeriod_s, refractoryPeriod_s],[0, 1])
area([0, refractoryPeriod_s],[1, 1], 'FaceAlpha', 0.1)



