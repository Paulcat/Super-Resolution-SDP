% Test Tprod2 (fft-based generalized Toeplitz product)

addpath('../../toolbox_signal')
addpath('../../toolbox_moment')

d = 2;


% 1- fft-based
% ------------
n = randi([8 12], [1,d]);
r = 2;
m = 2*n-1;

T = rand([m, 1]) + 1i*rand([m, 1]); % no need of Hermitianity for testing purpose
X = rand([n, r]) + 1i*rand([n, r]);

% Fast Toeplitz product
TX = Tprod2(T,X);


% Toeplitz helpers
[U1,V1] = meshgrid(1:n(1));
toepl1 = @(s) double(U1-V1==s);
[U2,V2] = meshgrid(1:n(2));
toepl2 = @(s) double(U2-V2==s);



% 2- full product
% ---------------
momo  = gencolex(n-1);
order = fft_to_colex(n-1);

Tmat = zeros(prod(n));
T    = T(order);
for i=1:size(momo,1)
    k    = momo(i,:);
    id   = k + n;
    Tmat = Tmat + T( id(1), id(2) ) * kron( toepl2(k(2)), toepl1(k(1)) );
end

TmatX = zeros([n, r]);
for i=1:r
    Xi           = X(:,:,i);
    TmatX(:,:,i) = reshape( Tmat*Xi(:), [n, 1] ); % for consistency with 1D
end



%
fprintf('\n\n\t---- Testing Tprod, and orderings ----\n');
fprintf('\t   Error (should be 0):%d\n\n', norm(TmatX(:) - TX(:), 'fro'));


