function train = poissonSpikeTrain(rate, rate_time, refractory)%, func)
% All inputs must be in units of seconds!

if (nargin < 2)
    error('poisson_spike_train: usage poisson_spike_train(rate, duration)')
end
if (nargin < 3)
    refractory = 0.001; % 1ms
end
%if (nargin < 4)
%    func = @(x)sigmf(x,[2/mean(rate),mean(rate)])*max(rate); % 1ms
%end

if (any(rate) < 0) % || any(rate(2:end).*diff(rate_time) > 1) )
    error('poisson_spike_train: rate must be a non-negative vector')
end

while any(rate(2:end).*diff(rate_time) > 1)
    %rate = repmat(rate, 2, 1); rate = rate(:)
    %M = ceil( max(rate(2:end).*diff(rate_time)) )  
    %rate = interp(rate, M);
    rate = interp1(1:length(rate), rate, 1:0.5:length(rate), 'linear');
    rate_time = interp1(1:length(rate_time), rate_time, 1:0.5:length(rate_time), 'linear');
end

if (~isscalar(duration) || duration < 0)
    error('poisson_spike_train: duration must be a non-negative scalar')
end

%%
%rate = func(rate);

%%
train_ = [];
for ii = 2:length(rate_time)
    interval = rate_time(ii)-rate_time(ii-1);
    nExpected = mean(rate((ii-1):ii))*interval;
    
    if rand(1) < nExpected
        train_ = [train_ rate_time(ii-1) + (0.25*randn(1)+0.5)*interval ];
    end
end

min_dt = 0;
while min_dt < refractory
    [min_dt min_idx] = min(diff(train_));
    train_(min_idx+1) = [];
end
train = train_;

%%
if 0
rate = (rate - min(rate))/(max(rate)-min(rate));
r = rand(size(rate));
train_ = rate_time(find(r>rate));


%dt = diff(train_);
min_dt = 0;
while min_dt < refractory
    [min_dt min_idx] = min(diff(train_));
    train_(min_idx+1) = [];
end
train = train_;
end
% if size(rate_time, 1) == 0
%     train = [];
% else
%     % Generate spike train <train_poiss> with poisson statistics
%     train = [];
%     t_next = -log(rand(1)) / randsample(rate, 1);
%     while t_next <= rate_time(end)
%         train = [train, t_next];
%         
%         delta = 0;
%         while delta < refractory
%             delta = - log(rand(1)) / rate(find(rate_time==t, 1, 'first'));
%         end
%         
%         t_next = t_next + delta;
%     end
% end
