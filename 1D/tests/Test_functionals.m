% Test functionals

addpath('../');
addpath('../toolbox_calculus/');
addpath('../../toolbox');
addpath('../../toolbox_signal/');
addpath('../../toolbox_moment/');



fc = 10;
n  = 2*fc + 1;
L  = 64;

lambda = 1e-3;
rho    = 1;
gam    = sqrt(1/L) * ( L < Inf ) + 1 * ( L == Inf );

model.problem = 'subsampled-convolution';
model.kernel = 'Gaussian';
model.kparam = 0.04;
model.grid   = 'lattice';
model.gsize  = L;
model.radius = 0.3;

[A,AS,w] = compute_A(fc,model);



y0 = rand( [L, 1] );
F = compute_fobj (fc,y0,lambda,rho,gam,A   );
G = compute_fgrad(fc,y0,lambda,rho,gam,A,AS);
G = @(U) 2 * G(U,U);


reshc = @(u,p) reshape(u(1:length(u)/2), [length(u)/(2*p),p]) + ... 
    1i*reshape(u(length(u)/2+1:end), [length(u)/(2*p),p]);
flatc = @(Z) [real(Z(:)); imag(Z(:))];

r  = 5;
U0 = rand( [n+1, r] ) + 1i * rand( [n+1, r] );

F = @(u) F(reshc(u,r));
G = @(u) flatc( G(reshc(u,r)) );
figure, checkgradient(F,G,flatc(U0));





