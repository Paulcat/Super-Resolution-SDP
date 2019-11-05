function [R,p] = perform_blasso_pen_sdp(fc,blasso,options)
%PERFORM_BLASSO_PEN_SDP solve "Toeplitz-penalized" BLASSO using SDP
%   Solve the following program:
%
%   min t + tr(R)/m + 1/(2*la)*||y-Az||^2 + 1/(2*rho)*||R - Toep(R)||^2
%   s.t.    |R  |z|
%           |___|_| >= 0
%           |z* |t|
%
%   Inputs: 
%       fc    - cutoff frequency (R is of size (2*fc+1))
%       blasso - blasso.y, etc... specifies the quantity composing the
%       objective

n = 2*fc + 1;
m = n;

solver = getoptions(options,'solver','primal');

y   = blasso.y;
la  = blasso.la;
rho = blasso.rho;
gam = blasso.ga;
C0  = blasso.f0;
A   = blasso.A;
As  = blasso.As;

l = length(y);

cvx_solver SDPT3 % SeDuMi | SDPT3
cvx_precision high


switch solver
    case 'primal'
        
        [X,Y] = meshgrid(1:m);
        toepl = @(s) double(X-Y==s);
        
        cvx_begin sdp %quiet
        variable R(m+1,m+1) hermitian
        %variable R1(m,m)    hermitian toeplitz
        variable z(n)       complex
        dual variable Q
        Q : R >= 0
        %
        T = zeros(m);
        for j=-m+1:m-1
            nj = abs(abs(j) - m); % number of diagonal elements
            dj = sum(diag(R(1:m,1:m),j)) / nj;
            T = T + dj * toepl(j);
        end
        R(1:n,m+1) == z
        if la==0
            error('TODO');
        else
            minimize( C0 * ( 1/2*(T(1,1) + R(end,end)) + ...
                1/2/la  * gam * norm(y - A(z), 'fro')^2 + ...
                1/2/rho * norm(R(1:m,1:m) - T,'fro')^2 ) );
        end
        cvx_end
        
        p = 1/la * (y-A(z));

    case 'dual'
        cvx_begin sdp
        variable Q(m+1,m+1) hermitian
        variable p(l) complex
        dual variable R
        R : Q >= 0
        Q(m+1,m+1) == 1;
        Q(1:n,m+1) == As(p);
        trace(Q) == 2
        for j = 1:m-1
            sum(diag(Q,j)) == Q(m+1-j,m+1);
        end
        if la==0
            warning('check factor gam');
            maximize( sqrt(gam)*real(p'*y) - rho/2 * norm(Q(1:m,1:m)-1/m*eye(m),'fro')^2 );
        else
            warning('check factor gam');
            maximize( sqrt(gam)*real(p'*y) - la/2*gam*norm(p,'fro')^2 - rho/2 * norm(Q(1:m,1:m)-1/m*eye(m),'fro')^2 );
        end
        cvx_end
end


end

