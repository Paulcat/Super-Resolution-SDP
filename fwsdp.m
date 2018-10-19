function [U] = fwsdp(fc,blasso,options)
%FWSDP Frank-Wolfe solver for SDP-BLASSO
%   FWSDP(fc,blasso,options) solves the following semidefinite program:
%
%   min f = 1/2*(t+u0) + 1/2/lambda*|y-A(z)|_H^2 + 1/2/rho*|R - Toep(R)|^2
%        |R  |z|
%   s.t. |___|_| >= 0
%        |z* |t|
%   where u0 = mean( diag(R) )
%   and y, lambda, rho, H are specified in *blasso*
%
%   
%   U = FWSDP(fc,blasso,options) returns the matrix U such that R = UU*.
%
%   The SDP is solved using a Frank-Wolfe approach, described in
%       *A Low-Rank Approach to Off-The-Grid Sparse Super-Resolution*
%       P. Catala, V. Duval and G. Peyre


% Setting options
display  = getoptions(options, 'display', 'on');
maxIter  = getoptions(options, 'maxIter', 20  );
tol      = getoptions(options, 'tol'    , 1e-5);

opt_bfgs = set_bfgs_options(options);
opt_lmo  = set_lmo_options (options);





n = 2*fc + 1;
m = n; % hierarchy order



F = blasso.obj;
g = blasso.grad;
G = @(U) 2 * g(U,U);


% Display infos
% -------------
if strcmp(display, 'on')
    fprintf( '\n\nCalling Frank-Wolfe solver\n' )
    
    fprintf( '------------------------------------------------------------\n')
    fprintf( [  'varsize    \t: %dx%d\n'    ...
                'lambda     \t: %.2e\n'      ...
                'rho        \t: %.2e\n'      ...
             ], prod(n)+1, prod(n)+1, blasso.lambda, blasso.rho       )
    fprintf( 'tolerance (dgap): %.0e\n', tol    )
    fprintf( 'bfgs tolerance\t: %.0e (%d)\n', opt_bfgs.progTol, opt_bfgs.MaxIter)
    fprintf( 'PI tolerance  \t: %.0e (%d)\n', opt_lmo.tol     , opt_lmo.maxIter )
    
    fprintf( '------------------------------------------------------------\n')
    
    fprintf('IT   OBJ\t  DGAP\t    PI (TIME)\t LS\t BFGS (TIME)\n')
    fprintf( '------------------------------------------------------------\n')
end
% -------------




% *** 0. Initialization ***
% -------------------------
tic;
U0     = zeros( prod(m) + 1, 1 );
%E0     = fobj(U0,y,fc,lambda,rho,gam,A);
E0     = F(U0);
U      = U0;
%U1     = reshape( U(1:end-1,:), [m, size(U,2)] );
E      = E0;
D0     = 2 * E0;
niter  = 0;
v0     = ones(prod(n)+1, 1) / sqrt(prod(n)+1); % inital vector for power iterations
n_PI    = [];
n_BFGS  = [];

sqinvOm = [ n(1) * ones(prod(n), 1); 1 ]; % n(1) is a HACK: howto when dimensions are not equal??





% *** 1. LMO step ***
% -------------------
% precompute gradient for lmo
%[~,dzU] = fgrad(U,h0,fc,y,lambda,rho,gam,A,AS);
%TU    = Tproj2( MatToTen(U(1:end-1,:), m) );
%GfU = @(h) sqinvOm .* fgrad_pre(U, sqinvOm .* h,dzU,TU,fc,y,lambda,rho,gam);
GfU = blasso.grad0U_handle(U);
GLMO = @(h) sqinvOm .* GfU( sqinvOm .* h );

time1 = toc;
[eVecm, eValm, infos] = perform_LMO_step(opt_lmo, GLMO, v0);
time1 = toc - time1;

if eValm > 0
    warning('Minimal eigenvalue is positive');
    eVecm = zeros(size(eVecm));
end
eVecm = sqrt(D0) * (sqinvOm .* eVecm);

% duality gap
% UGfU = U'     * fgrad(U,U,    fc,y,lambda,rho,gam,A,AS);
% eGfe = eVecm' * fgrad(U,eVecm,fc,y,lambda,rho,gam,A,AS);
UGfU = U'     * g(U, U    );
eGfe = eVecm' * g(U, eVecm);
gap    = real( trace(UGfU) - eGfe );




% Display infos
% -------------
if strcmp(display, 'on')
    fprintf('%-2d  %-+.4e  %-+.2e  -\t -\t -\n', niter, E(end), gap)
end

niter = 0;
while (gap >= tol && niter < maxIter)
    
    n_PI = [n_PI; infos.niter];
    
    % *** 2. Line-search ***
    % ----------------------
    [mu,nu] = perform_linesearch_step(fc,U,eVecm,blasso);
    
    
    % *** 3. FW update ***
    % --------------------
    U = [sqrt(mu)*U, sqrt(nu)*eVecm];
    %isToeplitz(fc,U(1:end-1,:))
    
    
    % *** 4. BFGS step ***
    % --------------------
    time2 = toc;
    [U,flag,output] = perform_bfgs_step(U,F,G,opt_bfgs);
    n_BFGS = [n_BFGS; output.iterations];
    time2 = toc-time2;
    
    %isToeplitz(fc,U(1:end-1,:))
    
    % One iteration done: update monitors
    niter = niter + 1;
    E     = [E; F(U)];
    
    % Display infos
    % -------------
    if strcmp(display,'on')
        fprintf('%-2d  %-+.4e\t %-+.2e  %-4i (%4.1f)\t %-.0e\t %-4i (%4.1f)\n', niter, E(end), gap, n_PI(end), time1, nu/mu, n_BFGS(end), time2)
    end
    % -------------
    
    
    % *** 1. LMO step ***
    % -------------------
    GfU = blasso.gradU_handle(U);
    
    time1 = toc;
    [eVecm, eValm, infos] = perform_LMO_step(opt_lmo, GfU, v0);
    time1 = toc - time1;
    
    if eValm > 0
        warning('Minimal eigenvalue is positive');
        eVecm = zeros(size(eVecm));
    end
    eVecm = sqrt(D0) * (sqinvOm .* eVecm);
    
    UGfU = U'     * g(U, U    );
    eGfe = eVecm' * g(U, eVecm);
    gap    = real( trace(UGfU) - eGfe );
    
end
time = toc;

info.E = [E; F(U)];
info.time = time;
info.nPI = sum(n_PI);
info.nBFGS = sum(n_BFGS);

if strcmp(display, 'on')
    fprintf('-------------------------------------------------------------\n')
    fprintf('DGAP: %+.2e\n', gap/E0)
    fprintf('OPTVAL: %+.5e\n', F(U))
    fprintf('TIME: %f\n', time)
end


end


function opt_bfgs = set_bfgs_options(options)
opt_bfgs.on              = 1;
opt_bfgs.display         = 'off';                                    % off | final | iter | full | excessive
opt_bfgs.optTol          = 1e-16;                                    % first-order optimality (exitflag 1)
opt_bfgs.progTol         = getoptions(options, 'bfgsProgTol', 1e-8); % parameters change (exitflag 2)
opt_bfgs.MaxFunEvals     = 100000;
opt_bfgs.MaxIter         = getoptions(options, 'bfgsMaxIter', 500);  % (exitflag 0)
opt_bfgs.Method          = 'lbfgs';
opt_bfgs.DerivativeCheck = 'off';
opt_bfgs.Corr            = 15;
opt_bfgs.Damped          = 0;
end

function opt_lmo = set_lmo_options(options)
opt_lmo.tol     = getoptions(options, 'lmoTpm', 1e-16);
opt_lmo.maxIter = getoptions(options, 'lmoMaxIter', 1000); 
end
