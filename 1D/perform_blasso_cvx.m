function [R,p,z] = perform_blasso_cvx(fc,blasso,options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n = 2*fc + 1;
m = n;

solver = getoptions(options, 'solver', 'primal');

[V,U] = meshgrid(1:m, 1:m);
toepl = @(s) double(U-V==s);

y   = blasso.y;
la  = blasso.lambda;
gam = blasso.ga;
C0  = blasso.normalization;
A   = blasso.A;


cvx_solver    mosek % SDPT3 | SeDuMi | mosek
cvx_precision high



%L = chol(inv(real(AS*A))); % real = HACK?

switch(solver)
    case 'primal'
        cvx_begin sdp %quiet
        variable R(m+1,m+1) hermitian;
        variable R1(m,m)    hermitian toeplitz;
        variable z(n)       complex;
        dual variable Q;
        Q : R >= 0
        %
        R(1:n,m+1) == z;
        R(1:m,1:m) == R1;
        if la==0
            error('TODO');
        else
            minimize( C0 * ( 1/2*(R1(1,1) + R(end,end)) + 1/2/la * gam * norm(y - A(z), 'fro')^2 ) );
        end
        cvx_end
        
        p = 1/la * (y - A(z)); % no AS*y, because PhiS is well implemented ??
        
    case 'cvx-dual' %TODO!!!!
        cvx_begin sdp %quiet
        variable Q(m+1,m+1) hermitian;
        variable Q1(m,m) hermitian;
        variable p(n) complex;
        dual variable R;
        R : Q >= 0;
        Q(m+1,m+1) == 1;
        Q(1:n,m+1) == p .* conj(w);
        Q(1:m,1:m) == Q1;
        Q(n+1:m,m+1) == 0;
        %         trace(Q) == 2;
        %         for j = 1:m-1,
        %             sum(diag(Q,j)) == Q(m+1-j,m+1);
        %         end
        for j=0:m-1
            T = topl(j);
            T(:)' * Q1(:) == double( j==0 );
        end
        if lambda==0
            maximize(real( p' * y ))
        else
            %maximize( real(p'*y) - lambda/2 * norm(p,'fro')^2 );
            %minimize( 1/2*norm(y/lambda-p,'fro')^2 );
            %minimize( lambda/2*norm(y/lambda-p,'fro')^2 - lambda/2*norm(y/lambda,'fro')^2);
            minimize( 1/2*sqrt(ga)*norm(L*A'*y/lambda - L*p)^2,'fro');
        end
        cvx_end
        
        p = A*p;
        


end

