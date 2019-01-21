% * ******************************* *
% * Super-resolution examples in 2D *
% * ******************************* *

clear all
path(pathdef)
addpath('2D/')
addpath('2D/toolbox')

addpath('toolbox')

% setup cvx -- change to your cvx location
addpath('~/Documents/MATLAB/cvx');
run cvx_setup.m

% setup bfgs
% reference: 
% M. Schmidt. minFunc: unconstrained differentiable multivariate 
% optimization in Matlab. 
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.
addpath(genpath([pwd '/minFunc/']));

%% synthetic data

s  = 6; % sparsity

x0 = [.4538 .8314; .4324 .8034; .4235 .3993; .1332 .5269; .3734 .4168; .3909 .6569];
a0 = [.7559; -.6160; .9690; .2681; .6657; -.7876];

% x0 = rand(s,2);
% a0 = -1 + 2*rand(s,1)

x0pos = x0(a0 >= 0, :); a0pos = a0(a0 >= 0); npos = length(a0pos);
x0neg = x0(a0 < 0, :); a0neg = a0(a0 < 0); nneg = length(a0neg);


% Fourier operator
fc     = [30, 30];
n      = 2*fc + 1;
[F,FS] = FourierOperator(fc);


%%  * Example 1: Gaussian convolution * 

sig = .02; % variance of isotropic Gaussian kernel
model.fop    = 'convolution';
model.kernel = 'Gaussian';
model.kparam = sig;
model.grid   = 'none'; % no subsampling grid

% forward operator
[A,AS] = approximationOperator(fc,model);
Phi  = @(x1,x2,a) A ( reshape( F(x1,x2) * a, [n,1] ) );
PhiS = @(dx,dy,p) FS( dx, dy, AS(p) );


% measurements
y0    = Phi(x0(:,1),x0(:,2),a0);
sigma = 1e-5;
wn    = randn(size(y0)) + 1i*randn(size(y0));
y     = y0 + sigma * norm(y0) * fft2(wn);

% Display
N  = 512;
dN = (0:N-1)'/N; % display grid step

clf, hold on
surf(dN,dN,FS(dN,dN,y),'LineStyle','none'), view(2);
contour3(dN,dN,FS(dN,dN,y0),10,'k','LineWidth',1.5);
stem3(x0pos(:,2),x0pos(:,1),ones(npos,1),'r.','MarkerSize',20);
stem3(x0neg(:,2),x0neg(:,1),ones(nneg,1),'c.','MarkerSize',20);
xlim([0.3,0.9]), ylim([0.05,0.6])
axis off

%% Solve BLASSO with FFW

% Blasso parameters
Cl     = norm( PhiS(dN,dN,y), 'inf'); % used for scaling lambda
lambda = Cl * 2e-3;
rho    = 1e3;
gam    = 1; % constant for Hilbert space norm (usual norm here since no grid)

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
p = 1/lambda * ( y - A( MatToTen( U(1:end-1,:) * U(end,:)', n) ) );

%% Support reconstruction

% Dual polynomial
eta = PhiS(dN,dN,p);

clf, hold on
surf(dN,dN,eta,'LineStyle','none');
view(2), axis tight, colorbar
stem3(x0pos(:,2),x0pos(:,1),ones(npos,1),'r.');
stem3(x0neg(:,2),x0neg(:,1),a0neg,'c.');


% support extraction procedure
U1 = U(1:end-1,:); % matrix main block
x1 = extract_spikes(U1,fc);

% amplitude recovery
[Yf,Xf] = FreqMesh(fc);
[~,~,Amat,ASmat] = approximationOperator(fc,model);
flat = @(x) x(:);
F_e = exp( -2i*pi*( Xf(:) * x1(:,1)' + Yf(:) * x1(:,2)' ) );
s0 = sign(real(F_e' * flat(AS(p)) ));
Phi_e = Amat * F_e;
PhiS_e = F_e' * ASmat;
a1 = real( Phi_e\y(:) - lambda * pinv(PhiS_e*Phi_e)*s0 );

x1pos = x1(a1 >= 0,:); a1pos = a1(a1 >= 0);
x1neg = x1(a1 < 0 ,:); a1neg = a1(a1 < 0 );

clf, hold on
stem3(x0pos(:,2),x0pos(:,1),a0pos,'ro','MarkerSize',18,'LineWidth',3 );
stem3(x0neg(:,2),x0neg(:,1),a0neg,'co','MarkerSize',18,'LineWidth',3 );
stem3(x1pos(:,2),x1pos(:,1),a1pos,'rx','MarkerSize',18,'LineWidth',3 );
stem3(x1neg(:,2),x1neg(:,1),a1neg,'cx','MarkerSize',18,'LineWidth',3 );



%% * Example 2: Subsampled Gaussian convolution *

sig = .02; % variance of isotropic Gaussian kernel
L   = [64,64]; % subsampling grid size
LX  = (0:L(1)-1)'/L(1);
LY  = (0:L(2)-1)'/L(2);

model.fop    = 'subsampled-convolution';
model.kernel = 'Gaussian';
model.kparam = sig;
model.grid   = 'lattice';
model.gsize  = L;

% forward operator
[A,AS] = approximationOperator(fc,model);
Phi    = @(x1,x2,a) A( reshape( F(x1,x2)*a, [n,1] ) );
PhiS   = @(dx,dy,p) FS( dx, dy, AS(p) );

% measurements
y0    = Phi(x0(:,1),x0(:,2),a0);
sigma = 1e-2;
wn    = randn(size(y0));
y     = y0 + sigma * norm(y0) * wn;

clf, hold on
surf(LX,LY,real(y),'LineStyle','none'), view(2);
%contour3(LX,LY,real(y0),9,'k','LineWidth',1.5);
stem3(x0pos(:,2),x0pos(:,1),2*ones(npos,1),'r.','MarkerSize',20);
stem3(x0neg(:,2),x0neg(:,1),ones(nneg,1),'c.','MarkerSize',20);
xlim([0.3,0.9]), ylim([0.05,0.6])
axis off

%% Solve BLASSO with FFW

% Blasso parameters
Cl     = norm( PhiS(dN,dN,y), 'inf'); % used for scaling lambda
lambda = Cl * 1e-3;
rho    = 1e4;
gam    = sqrt(1/prod(L)); % constant for Hilbert space norm

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
p = 1/lambda * ( y - A( MatToTen( U(1:end-1,:) * U(end,:)', n) ) );

%% Support reconstruction

% Dual polynomial
eta = PhiS(dN,dN,p);

clf, hold on
surf(dN,dN,eta,'LineStyle','none');
view(2), axis tight, colorbar
stem3(x0pos(:,2),x0pos(:,1),ones(npos,1),'r.');
stem3(x0neg(:,2),x0neg(:,1),a0neg,'c.');


% support extraction procedure
U1 = U(1:end-1,:); % matrix main block
x1 = extract_spikes(U1,fc);

% amplitude recovery
[Yf,Xf] = FreqMesh(fc);
[~,~,Amat,ASmat] = approximationOperator(fc,model);
flat = @(x) x(:);
F_e = exp( -2i*pi*( Xf(:) * x1(:,1)' + Yf(:) * x1(:,2)' ) );
s0 = sign(real(F_e' * flat(AS(p)) ));
Phi_e = Amat * F_e;
PhiS_e = F_e' * ASmat;
a1 = real( Phi_e\y(:) - lambda * pinv(PhiS_e*Phi_e)*s0 );

x1pos = x1(a1 >= 0,:); a1pos = a1(a1 >= 0);
x1neg = x1(a1 < 0 ,:); a1neg = a1(a1 < 0 );

clf, hold on
stem3(x0pos(:,2),x0pos(:,1),a0pos,'ro','MarkerSize',18,'LineWidth',3 );
stem3(x0neg(:,2),x0neg(:,1),a0neg,'co','MarkerSize',18,'LineWidth',3 );
stem3(x1pos(:,2),x1pos(:,1),a1pos,'rx','MarkerSize',18,'LineWidth',3 );
stem3(x1neg(:,2),x1neg(:,1),a1neg,'cx','MarkerSize',18,'LineWidth',3 );