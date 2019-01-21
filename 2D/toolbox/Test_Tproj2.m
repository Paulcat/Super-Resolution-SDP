% Test Tproj2 (fft-based generalized Toeplitz projection)

addpath('../../toolbox_signal')
addpath('../../toolbox_moment')

d = 2;

n = 5*ones(1,d);
r = 2;

[X,Y] = meshgrid(1:n(1));
toepl = @(s) double(X-Y==s);


% 1- fft-based projection
% -----------------------
U = rand([n, r]);
T = Tproj2(U);

momo  = gencolex(n-1);
order = fft_to_colex(n-1);

T     = T(order);
fastT = zeros(prod(n));
for i=1:size(momo,1)
    k  = momo(i,:);
    id = k+n;
    fastT = fastT + T( id(1), id(2) ) * kron( toepl(k(2)), toepl(k(1)) );
%     imagesc( kron( toepl(k(2)), toepl(k(1)) ) ); drawnow;
%     imagesc( kron( toepl(k(2)), toepl(k(1)) ) ); drawnow;
%     imagesc( kron( toepl(k(2)), toepl(k(1)) ) ); drawnow;
%     imagesc( kron( toepl(k(2)), toepl(k(1)) ) ); drawnow;
%     imagesc( kron( toepl(k(2)), toepl(k(1)) ) ); drawnow;
end
T     = T(order);


% 2- naive projection
% -------------------
U2 = reshape(U, [prod(n), r]);
M  = U2 * U2';

momo = gencolex(n-1);

naiveT = zeros(prod(n));
for i=1:size(momo,1)
    k = momo(i,:);
    Theta = kron( toepl(k(2)), toepl(k(1)) );
    Theta_n = Theta / norm(Theta,'fro');
    naiveT = naiveT + trace(Theta_n'*M) * Theta_n;
end

fprintf('\n\n\t---- Testing Tproj, and orderings ----\n');
fprintf('\t  Error (should be 0): %d\n\n', norm(fastT-naiveT,'fro') / norm(naiveT,'fro'));
