function [U,info] = ffw(fc,blasso,options)
%FFW Frank-Wolfe solver for SDP-BLASSO
%   FWSDP(fc,blasso,options) solves the following semidefinite program:
%
%   min f = 1/2*(t+tr(R)/m) + 1/2/la*|y-A(z)|_H^2 + 1/2/rho*|R - Toep(R)|^2
%        |R  |z|
%   s.t. |___|_| >= 0
%        |z* |t|
%
%   where y, lambda, rho, H are specified in *blasso*
%
%   
%   U = FWSDP(fc,blasso,options) returns the matrix U such that R = UU*.
%
%   The SDP is solved using a Frank-Wolfe approach, described in
%       *A Low-Rank Approach to Off-The-Grid Sparse Super-Resolution*
%       P. Catala, V. Duval and G. Peyre


debug = getoptions(options,'mode_debug',0);


% Setting options
display  = getoptions(options, 'display', 'on');
maxIter  = getoptions(options, 'maxIter', 20  );
tol      = getoptions(options, 'tol'    , 1e-5);

opt_bfgs = set_bfgs_options(options);
opt_lmo  = set_lmo_options (options);


% count number of fft calls
global nfft
nfft = 0;
 
n = 2*fc + 1;
m = n; % hierarchy order



F  = blasso.obj  ;
g  = blasso.grad ;
G  = @(U) 2 * g (U,U);





% Display infos
% -------------
if strcmp(display, 'on')
    fprintf( '\n\nCalling Frank-Wolfe solver\n' )
    
    fprintf( '------------------------------------------------------------\n')
    fprintf( [  'varsize    \t: %dx%d\n'    ...
                'lambda     \t: %.2e\n'      ...
                'rho        \t: %.2e\n'      ...
             ], prod(n)+1, prod(n)+1, blasso.la, blasso.rho       )
    fprintf( 'tolerance (dgap): %.0e\n', tol    )
    fprintf( 'bfgs tolerance\t: %.0e (%d)\n', opt_bfgs.progTol, opt_bfgs.MaxIter)
    fprintf( 'PI tolerance  \t: %.0e (%d)\n', opt_lmo.tol     , opt_lmo.maxIter )
    
    fprintf( '------------------------------------------------------------\n')
    
    fprintf('IT   OBJ\t  DGAP\t    PI   (TIME)\t LS\t BFGS (TIME)\n')
    fprintf( '------------------------------------------------------------\n')
end
% -------------






% * Initialization *
% -------------------------
tic;
U0     = zeros( prod(m) + 1, 1 ); U = U0; % primal variable
%mu0L   = zeros( prod(m), 1); muL = mu0L;
%mu0T   = zeros( 2*prod(m)-1, 1); muT = mu0T; % dual variables
E0     = F(U0); E = E0;
D0     = 2 * E0;
%D0 = 2 * prod(m) * E0;
niter  = 0;
v0     = ones(prod(n)+1, 1) / sqrt(prod(n)+1); % initial vector for PI
n_PI   = [];
n_BFGS = [];

Om = [ n(1) * ones(prod(n), 1); 1 ]; % TODO: handle non-equal dimensions
%sqinvOm = ones(prod(n)+1,1);

% % Display infos
% % -------------
% if strcmp(display, 'on')
%     fprintf('%-2d  %-+.4e  %-+.2e  -\t -\t -\t -\n', niter, E(end), gap)
% end

Gval = [];
gap = Inf;


%clf;
while (gap >= tol && niter < maxIter)
    
    % * 1. LMO *
    Glmo = blasso.gradHandle(U);
    GLMO = @(h) Om .* Glmo(Om .* h);
    
    % added for the CGAL algorithm
    % NOTE useless at 1st step since dual variable is initialized at
    % zero... (?)
    % ----------------------------
    %GLU = blasso.LgradU_handle(muL,muT,U);
    %GfU2 = @(h) GfU(h) + GLU(h);
    %GLMO2 = @(h) sqinvOm .* GfU2( sqinvOm .* h );
    % ----------------------------
    
    time1 = toc;
    [eVecm, eValm, eVecM, infos] = perform_LMO_step(opt_lmo, GLMO, v0);
    time1 = toc - time1;
    
    Gval = [Gval, real(eVecM'*g(U,eVecM))];
    %GLMO2 = @(h) GLMO(h)/Gval(end);
    %eVecm = perform_LMO_step(opt_lmo,GLMO2, v0);
    n_PI = [n_PI; infos.niter];
    
    if Gval(end)/norm(U,'fro') < 1e-10
        fprintf('gradient very close to zero: %d\n', Gval(end));
        break;
    end
    
    if debug
        dif = @(a,b) norm(a-b,'fro')/norm(b,'fro');
        
        D = test_grad(U,eVecm,blasso);
        DOm = diag(Om) * D * diag(Om);
        [VV,DD] = eig(DOm);
        
        % call perform_LMO with explicit gradient multiplication
        GLMO2 = @(h) DOm*h;
        [evec2,eval2] = perform_LMO_step(opt_lmo,GLMO2,v0);
        
        nn = max(imag(diag(DD)));
        if nn > 1e-13
            error('gradient not Hermitian, max imaginary part of eigenvalues: %d\n',nn);
        end
        
        [meval,id] = min(real(diag(DD)));
        mevec = VV(:,id);
        mm = norm(DOm*mevec - meval*mevec,'fro');
        if mm > 1e-12
            error('problem with gradient eigen decomposition, |(A-xI)v| = %d\n', mm);
        end
        
        fprintf('\t checking PI\n')
        fprintf('\t\t- fast vs matlab eig: %d\n', abs(eValm-meval)/abs(meval));
        fprintf('\t\t- calling PI with explicit multiplication: evec: %d, eval: %d\n', dif(eVecm,evec2), dif(eValm,eval2));
        fprintf('\t\t- verify eigen relation: %d\n', norm(DOm*eVecm - meval*eVecm,'fro')/abs(meval)); 
        %fprintf('eValm: %.5d, trace: %.5d\n', eValm, -1/2 * real(trace(U'*G(U))) );
    end
    
    if eValm > 0
        fprintf('\n Stopping: gradient is positive, minimal eigenvalue = %d\n', eValm);
        break;
    end
    eVecm = sqrt(D0) * (Om .* eVecm);
    
    % compute duality gap
    UGU = U'     * g(U, U    );
    eGe = eVecm' * g(U, eVecm);
    gap  = real( trace(UGU) - eGe );
    
    
    % * 2. FW update (with optional linesearch) *
    [mu,nu] = perform_linesearch_step(fc,U,eVecm,blasso.ls);
    if mu == 0
        U = sqrt(nu)*eVecm;
    elseif nu == 0
        U = sqrt(mu)*U;
    else
        U = [sqrt(mu)*U, sqrt(nu)*eVecm];
    end 
    E = [E;F(U)];
    
    % TODO: est-ce que si à une étape de l'algo on satisfait la contrainte
    % toeplitz, ca reste vrai jusqu'à la fin de l'algo??
    
    
    % * 3. BFGS step *
    time2 = toc;
    [U,flag,output] = perform_bfgs_step(U,F,G,opt_bfgs);
    n_BFGS = [n_BFGS; output.iterations];
    time2 = toc-time2;
    
    
    % added for CGAL algorithm
    % ------------------------
    % *** 5. Dual update ***
    %theta = 1/(niter+1);
    %muL = [muL, sqrt(theta) * U(1:end-1,:)];
    %muT = muT + theta * Tproj1(U(1:end-1,:));
    % ------------------------
    

    
    % One iteration done: update monitors
    niter = niter + 1;
    E = [E;F(U)];
    
    % Display infos
    % -------------
    if strcmp(display,'on')
        fprintf('%-2d  %-+.4e\t %-+.2e  %-4i (%4.1f)\t %-.0e\t %-4i (%4.1f)\n', niter, E(end), gap, n_PI(end), time1, nu/mu, n_BFGS(end), time2)
    end
    % -------------
    
    E2 = E(1:2:end);
    if abs(E2(end) - E2(end-1))/abs(E2(end)) < 1e-10
        % stopping criterion in terms of objective decrease
        sprintf('\n Stopping: energy did not decrease at last step\n')
        break;
    end
    
end
time = toc;

if niter == maxIter-1
    fprintf('\n Stopping: maximum number of iterations reached\n');
elseif gap < tol
    fprintf('\n Stopping: criterion is satisfied\n');
end


info.E = [E; F(U)];
info.G = Gval;
info.time = time;
info.nPI = sum(n_PI);
info.nBFGS = sum(n_BFGS);
info.nfft  = nfft;

if debug
    clf,
    subplot(1,2,1), semilogy(1:2*(niter+1), log10(info.E)); title('E');
    subplot(1,2,2), semilogy(1:niter+1, log10(info.G)); title('G');
end


%U(:,1) = []; % first columns is only zeros


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
opt_bfgs.numDiff         = 0; % use-provided gradient
end

function opt_lmo = set_lmo_options(options)
opt_lmo.tol     = getoptions(options, 'lmoTol', 1e-16);
opt_lmo.maxIter = getoptions(options, 'lmoMaxIter', 1000); 
end
