% * Super-resolution examples in 2D *
% * ******************************* *

clear all
path(pathdef)
addpath('../')
addpath('../2D/')
addpath('../2D/toolbox')

addpath('../toolbox')

% setup cvx -- change to your cvx location
addpath('~/Documents/MATLAB/cvx');
run cvx_setup.m

% setup bfgs
% reference: 
% M. Schmidt. minFunc: unconstrained differentiable multivariate 
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

%% Synthetic data

s  = 6; % sparsity

x0 = [.4538 .8314; .4324 .8034; .4235 .3993; .1332 .5269];%; .3734 .4168; .3909 .6569];
a0 = [.7559; -.6160; .9690; .2681];%; .6657; -.7876];

x0pos = x0(a0 >= 0, :); a0pos = a0(a0 >= 0); npos = length(a0pos);
x0neg = x0(a0 < 0, :); a0neg = a0(a0 < 0); nneg = length(a0neg);

% Fourier operator
fc     = [30, 30];
n      = 2*fc + 1;
[F,Fs] = FourierOperator(fc);

% display grid
N  = 512;
dN = (0:N-1)'/N; % display grid step

%% * Example 1: Ideal low-pass filtering

clear model
model.fop = 'convolution';

model.kernel.type = 'Dirichlet';

gam = 1; % constant for Hilbert space norm

%% * Example 2: Gaussian convolution * 

clear model
model.fop = 'convolution';

% kernel (isotropic Gaussian)
sig = .04;
model.kernel.type = 'Gaussian';
model.kernel.cov = sig;

gam = 1; % constant for Hilbert space norm

%% * Example 3: Subsampled Gaussian convolution

clear model
model.fop = 's-convolution';

% kernel (isotropic Gaussian)
sig = .02;
model.kernel.type = 'Gaussian';
model.kernel.cov = sig;

% subsampling grid
L   = [64,64];
LX  = (0:L(1)-1)'/L(1);
LY  = (0:L(2)-1)'/L(2);
model.grid.shape   = 'lattice';
model.grid.size  = L;
gam = sqrt(1/prod(L)); % constant for Hilbert space norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM

%% Measurements

% forward operator
[A,As,AsA] = approximationOperator(fc,model);
%Phi  = @(x1,x2,a) A ( reshape( F(x1,x2) * a, [n,1] ) );
Phi = @(x1,x2,a) A( F(x1,x2)*a );
PhiS = @(dx,dy,p) Fs( dx, dy, As(p) );

% measurements
y0    = Phi(x0(:,1),x0(:,2),a0); % noiseless observations
sigma = 0;
wn    = randn(size(y0)) + 1i*randn(size(y0));
y     = y0 + sigma * norm(y0) * fft2(wn);

figure(1), clf, hold on
if strcmp(model.fop, 'convolution')
    % for convolution, y lives in Fourier space: plot Fs(y)
    surf(dN,dN,Fs(dN,dN,y),'LineStyle','none'), view(2);
    contour3(dN,dN,Fs(dN,dN,y0),10,'k','LineWidth',1.5);
else
    surf(LX,LY,reshape(real(y),L), 'linestyle','none');
end
stem3(x0pos(:,2),x0pos(:,1),a0pos,'r.','MarkerSize',20);
stem3(x0neg(:,2),x0neg(:,1),abs(a0neg),'c.','MarkerSize',20);
xlim([0.3,0.9]), ylim([0.05,0.6])
view(2), colorbar, axis off

%% Solve BLASSO with FFW

% Blasso parameters

% scalings
Cl = norm( PhiS(dN,dN,y), 'inf'); % scaling for lambda
Cr = 1; % scaling for rho

la = Cl * 1e-3; % sparse regularization
rho    = Cr * 1e2; % toeplitz penalization

% main parameters
blasso.y      = y;
blasso.la = la;
blasso.rho    = rho;
blasso.ga     = gam;
blasso.A      = A;
blasso.As     = As;

% objective and gradient
blasso.obj        =      fobj (   fc,y,la,rho,gam,As,AsA);
blasso.grad       =      fgrad(   fc,y,la,rho,gam,As,AsA);
blasso.gradHandle = @(U) fgradU(U,fc,y,la,rho,gam,As,AsA);
blasso.ls         = ls_coeffs(fc,y,la,rho,gam,As,AsA);

% FFW options
options.maxIter     = 10;
options.bfgsProgTol = 1e-10;
options.bfgsMaxIter = 500;
options.lmoTol      = 1e-10;
options.lmoMaxIter  = 1000;
options.tol         = 1e-6; % tolerance on dual gap

U = ffw(fc,blasso,options);
%p = 1/lambda * ( y - A( MatToTen( U(1:end-1,:) * U(end,:)', n) ) );
z = U(1:end-1,:)*U(end,:)';
p = 1/la * ( y - A(z) );

%% Support reconstruction

% Dual polynomial
eta = PhiS(dN,dN,p);

figure(1), clf, hold on
surf(dN,dN,eta,'LineStyle','none');
view(2), axis tight, colorbar
stem3(x0pos(:,2),x0pos(:,1),ones(npos,1),'r.');
stem3(x0neg(:,2),x0neg(:,1),a0neg,'c.');


% support extraction procedure
U1 = U(1:end-1,:); % matrix main block
%x1 = extract_spikes(U1,fc);
options.mode_debug = 0;
options.tol = 1e-3;
x1 = mvprony(U1*U1',2,options);

% amplitude recovery
[Yf,Xf] = FreqMesh(fc);
[~,~,~,Amat,ASmat] = approximationOperator(fc,model);
flat = @(x) x(:);
F_e = exp( -2i*pi*( Xf(:) * x1(:,1)' + Yf(:) * x1(:,2)' ) );
s0 = sign(real(F_e' * flat(As(p)) ));
Phi_e = Amat * F_e;
PhiS_e = F_e' * ASmat;
a1 = real( Phi_e\y(:) - la * pinv(PhiS_e*Phi_e)*s0 );

x1pos = x1(a1 >= 0,:); a1pos = a1(a1 >= 0);
x1neg = x1(a1 < 0 ,:); a1neg = a1(a1 < 0 );

figure(2), clf, hold on
stem3(x0pos(:,2),x0pos(:,1),a0pos,'ro','MarkerSize',18,'LineWidth',3 );
stem3(x0neg(:,2),x0neg(:,1),a0neg,'co','MarkerSize',18,'LineWidth',3 );
stem3(x1pos(:,2),x1pos(:,1),a1pos,'rx','MarkerSize',18,'LineWidth',3 );
stem3(x1neg(:,2),x1neg(:,1),a1neg,'cx','MarkerSize',18,'LineWidth',3 );

