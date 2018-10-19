% Test Tproj1

addpath('../../toolbox_signal')
addpath('../../toolbox_moment')

d = 1;

n = 5;
r = 2;

[X,Y] = meshgrid(1:n);
toepl = @(s) double(X-Y==s);


% 1- fft-base projection
% ----------------------
U = rand([n, r]) + 1i * rand([n, r]);
T = Tproj1(U);

momo  = gencolex(n-1);
order = fft_to_colex(n-1);

T     = T(order);
fastT = zeros(prod(n));
for i=1:size(momo,1)
    k  = momo(i,:);
    id = k+n;
    
    fastT = fastT + T(id) * toepl(k);
end
T = T(order);


% 2- naive projection
% -------------------
M = U*U';

momo = gencolex(n-1);

naiveT = zeros(n);
for i=1:size(momo,1)
    k       = momo(i,:);
    Theta   = toepl(k);
    Theta_n = Theta / norm(Theta,'fro');
    
    naiveT  = naiveT + trace(Theta_n'*M) * Theta_n;
end