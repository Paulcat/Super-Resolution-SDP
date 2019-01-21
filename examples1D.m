% * ******************************* *
% * Super-resolution examples in 1D *
% * ******************************* *

clear all
path(pathdef)
addpath('1D/')
addpath('1D/toolbox')

addpath('toolbox')

% setup cvx -- change to your cvx location
addpath('~/Documents/MATLAB/cvx')
run cvx_setup.m

% setup bfgs
% see: M. Schmidt. minFunc: unconstrained differentiable multivariate 
% optimization in Matlab. 
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.
addpath(genpath([pwd '/minFunc/']));

%% synthetic data

s = 4; % sparsity

x0 = [.4321; .6003; .3001; .6499];
a0 = [.2; -.5; .4; .6];

x0pos = x0(a0 >= 0, :); a0pos = a0(a0 >= 0);
x0neg = x0(a0 < 0,  :); a0neg = a0(a0 < 0);

% Fourier operator
fc = 30;
n = 2*fc + 1;
[F,FS] = FourierOperator(fc);


%% * Example 1: Foveation *

sig = @(x) .002 + .07*abs(x-.5);
L = 256;
dL = (0:L-1)'/L;

model.fop    = 'foveation';
model.kernel = 'Gaussian';
model.kparam = sig;
model.grid   = 'lattice';
model.gsize   = L;

% forward operator
[A,AS] = approximationOperator(fc,model);
Phi  = @(x,a)  A ( reshape( F(x) * a, [n,1] ) );
PhiS = @(dx,p) FS( dx, AS(p) );

% measurements
y0    = Phi(x0,a0);
sigma = 1e-3;
wn    = randn(size(y0));
y     = y0 + sigma * norm(y0) * wn;

clf, hold on
plot(dL,real(y),'LineWidth',2);
stem(x0pos,a0pos,'r.','MarkerSize',18);
stem(x0neg,a0neg,'b.','MarkerSize',18);

%% Solve BLASSO with FFW

% Blasso parameters
N  = 512;
dN = (0:N-1)'/N; % display grid step
Cl = norm( PhiS(dN,y), 'inf' );
lambda = Cl * 1e-3;
rho = 1e2;
gam = sqrt(1/prod(L));

[blasso.obj, blasso.f0] =      fobj (   fc,y,lambda,rho,gam,A   );
blasso.grad             =      fgrad(   fc,y,lambda,rho,gam,A,AS);
blasso.gradU_handle     = @(U) fgradU(U,fc,y,lambda,rho,gam,A,AS);
blasso.y                = y;
blasso.lambda           = lambda;
blasso.rho              = rho;
blasso.ga               = gam;
blasso.A                = A;
blasso.AS               = AS;

% FFW options
options.maxIter     = 10;
options.bfgsProgTol = 1e-9;
options.bfgsMaxIter = 500;
options.lmoTol      = 1e-10;
options.lmoMaxIter  = 1000;
options.tol         = 1e-6; % tolerance on dual gap

U = ffw(fc,blasso,options);
p = 1/lambda * ( y - A( U(1:end-1,:) * U(end,:)' ) );

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
x1 = extract_spikes(U1,fc);

% amplitude recovery
Xf = FreqMesh(fc);
[~,~,Amat,ASmat] = approximationOperator(fc,model);
flat = @(x) x(:);
F_e = exp( -2i*pi*( Xf(:) * x1(:)' ) );
s0 = sign(real(F_e' * flat(AS(p)) ));
Phi_e = Amat * F_e;
PhiS_e = F_e' * ASmat;
a1 = real( Phi_e\y(:) - lambda * pinv(PhiS_e*Phi_e)*s0 );

x1pos = x1(a1 >= 0,:); a1pos = a1(a1 >= 0);
x1neg = x1(a1 < 0 ,:); a1neg = a1(a1 < 0 );

clf, hold on
stem(x0pos,a0pos,'ro','MarkerSize',10,'LineWidth',2 );
stem(x0neg,a0neg,'co','MarkerSize',10,'LineWidth',2 );
stem(x1pos,a1pos,'rx','MarkerSize',15,'LineWidth',2 );
stem(x1neg,a1neg,'cx','MarkerSize',15,'LineWidth',2 );
xlim([0 1])
