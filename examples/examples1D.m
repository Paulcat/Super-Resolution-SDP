% * Super-resolution examples in 1D *
% * ******************************* *

clear all
path(pathdef)
addpath('../')
addpath('../1D/')
addpath('../1D/toolbox')

addpath('../toolbox')

% setup cvx -- change to your cvx location
addpath('~/Documents/MATLAB/cvx')
run cvx_setup.m

% setup bfgs
% see: M. Schmidt. minFunc: unconstrained differentiable multivariate 
% optimization in Matlab. 
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.
addpath('../minFunc/')
addpath('../minFunc/minFunc/')
addpath('../minFunc/minFunc/compiled')
addpath('../minFunc/minFunc/mex')
addpath('../minFunc/autoDif')
%addpath(genpath([pwd 'minFunc/']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP

%% synthetic data

s = 4; % sparsity

x0 = [.4321; .6003; .3001; .6499];
a0 = [.2; -.5; .4; .6];

x0pos = x0(a0 >= 0, :); a0pos = a0(a0 >= 0);
x0neg = x0(a0 < 0,  :); a0neg = a0(a0 < 0);

% Fourier operator
fc = 20;
n = 2*fc + 1;
[F,Fs] = FourierOperator(fc);


% display grid
N  = 512;
dN = (0:N-1)'/N;

%% * Example 1: Ideal low-pass filtering * 

clear model
model.fop = 'convolution';
model.kernel.type = 'Dirichlet';

% parameter for the Hilbert norm (in BLASSO)
gam = 1;

%% * Example 2: Gaussian convolution *

clear model
model.fop = 'convolution';
model.kernel.type = 'Gaussian';
model.kernel.cov = .02;

gam = 1;

%% * Example 3: Subsampled Gaussian convolution *

clear model
model.fop = 's-convolution';
model.kernel.type = 'Gaussian';
model.kernel.cov = .02;

% sampling grid
L = 64;
dL = (0:L-1)'/L;
model.grid.shape = 'lattice';
model.grid.size = L;

%% * Example 4: Foveation *

clear model
model.fop    = 'foveation';

% kernel (Gaussian with varying std)
model.kernel.type = 'Gaussian';
sig = @(x) .002 + .07*abs(x-.5);
model.kernel.cov = sig;

% sampling grid
L = 256;
dL = (0:L-1)'/L;
model.grid.shape  = 'lattice';
model.grid.size  = L;

% parameter for the Hilbert norm (in BLASSO)
gam = sqrt(1/prod(L));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM

%% Measurements

% forward operator
[A,As,AsA] = approximationOperator(fc,model);
Phi  = @(x,a)  A ( reshape( F(x) * a, [n,1] ) );
PhiS = @(dx,p) Fs( dx, As(p) );

% measurements
y0    = Phi(x0,a0);
sigma = 5e-4;
wn    = randn(size(y0));
y     = y0 + sigma * norm(y0) * wn;

clf, hold on

stem(x0pos,a0pos,'r.','MarkerSize',18);
stem(x0neg,a0neg,'b.','MarkerSize',18);

if strcmp(model.fop, 'convolution')
    % for convolution, y lives in Fourier space: plot Fs(y)
    plot(dN,Fs(dN,y),'LineWidth',2);
else
    plot(dL,real(y),'LineWidth',2);
end

%% Solve BLASSO with FFW

% Blasso parameters
Cl = norm( PhiS(dN,y), 'inf' );
Cr = 1;

la = Cl * 5e-3; % sparse regularization
rho = Cr * 1e1; % toeplitz penalization

% main parameters
blasso.As  = As;
blasso.AsA = AsA;
blasso.y   = y;
blasso.la  = la;
blasso.rho = rho;
blasso.ga  = gam;

% objective and gradient
blasso.obj        =      fobj (   n,y,la,rho,gam,As,AsA);
blasso.grad       =      fgrad(   n,y,la,rho,gam,As,AsA);
blasso.gradHandle = @(U) fgradU(U,n,y,la,rho,gam,As,AsA);
blasso.ls         = ls_coeffs(n,y,la,rho,gam,As,AsA);

% FFW options
options.maxIter     = 10;
options.bfgsProgTol = 1e-9;
options.bfgsMaxIter = 500;
options.lmoTol      = 1e-10;
options.lmoMaxIter  = 1000;
options.tol         = 1e-6; % tolerance on dual gap

U = ffw(fc,blasso,options);
p = 1/la * ( y - A( U(1:end-1,:) * U(end,:)' ) );

%% Support reconstruction

% Dual polynomial
eta = PhiS(dN,p);

clf, hold on
plot(dN,eta,'LineWidth',2);
stem(x0pos,a0pos,'r.');
stem(x0neg,a0neg,'b.');
axis tight

% support extraction procedure
U1 = U(1:end-1,:); % matrix main block
%x1 = extract_spikes(U1,fc);
options.mode_debug = 0;
options.tol = 1e-3;
x1 = mvprony(U1*U1',1,options);

% amplitude recovery
Xf = FreqMesh(fc);
[~,~,~,Amat,Asmat] = approximationOperator(fc,model);
flat = @(x) x(:);
F_e = exp( -2i*pi*( Xf(:) * x1(:)' ) );
s0 = sign(real(F_e' * flat(As(p)) ));
Phi_e = Amat * F_e;
PhiS_e = F_e' * Asmat;
a1 = real( Phi_e\y(:) - la * pinv(PhiS_e*Phi_e)*s0 );

x1pos = x1(a1 >= 0,:); a1pos = a1(a1 >= 0);
x1neg = x1(a1 < 0 ,:); a1neg = a1(a1 < 0 );

clf, hold on
stem(x0pos,a0pos,'ro','MarkerSize',10,'LineWidth',2 );
stem(x0neg,a0neg,'co','MarkerSize',10,'LineWidth',2 );
stem(x1pos,a1pos,'rx','MarkerSize',15,'LineWidth',2 );
stem(x1neg,a1neg,'cx','MarkerSize',15,'LineWidth',2 );
xlim([0 1])
