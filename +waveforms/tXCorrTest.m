
%% ---
dt = 0.05;
t = 0:dt:5;
D = 7
D2 = 22

z0 = zeros(100, 3);
z_shift = zeros(D, 3);
z_shift2 = zeros(D2, 3);
x_ = [sin(t); 0.5*sin(t); 0.5*sin(t+0.1)]';

x1 = [z; x_; z_shift; z_shift2; z];
x2 = [z; z_shift; x_; z_shift2; z];
x3 = [z; z_shift2; x_; z_shift; z];

% ---
T = [];
T(:,:,1) = x1;
T(:,:,2) = x2;
T(:,:,3) = x3;
[xc A S] = waveforms.tXCorr(T);

%% ---
figure; hold on;
plot(reshape(x1, 1, []), 'b');
plot(reshape(x2, 1, []), 'r');
plot(reshape(x3, 1, []), 'm');

n2 = S(2,1);
x2_shifted = [x2((n2+1):end, :); zeros(n2,size(x2, 2))];
plot(reshape( x2_shifted, 1, [])+0.01, 'g');

n3 = S(3,1);
x3_shifted = [x3((n3+1):end, :); zeros(n3,size(x3, 2))];
plot(reshape( x3_shifted, 1, [])+0.02, 'c');


