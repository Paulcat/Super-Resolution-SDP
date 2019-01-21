function [R,p] = perform_blasso_sdp(fc,blasso,options)
%PERFORM_BLASSO_CVX Solve the Blasso using cvx-interfaced interior points
%   [R,p] = PERFORM_BLASSO_CVX(fc,blasso,options) returns the moment matrix
%   R and the dual vector p of a solution of the BLASSO.
%
%   options.solver: either 'primal' or 'dual'
%
%   The idea of solving the problem using a SDP lifting is introduced in:
%       *Towards a Mathematical Theory of Super-Resolution*
%       E. J. Candes and C. Fernandez-Granda
%   See also
%       *Positive Trigonometric Polynomials and Sig. Proc. Applications*
%       B. Dumitrescu
%
%   Copyright (c) 2013 Gabriel Peyre

n = 2*fc + 1;
m = n;

solver = getoptions(options, 'solver', 'cvx-primal');

[V,U] = meshgrid(1:m, 1:m);
toepl = @(s) double(U-V==s);

indicator = Pad2(ones(n),m);
I         = indicator(:) > 0;

y   = blasso.y;
la  = blasso.lambda;
gam = blasso.ga;
C0  = blasso.normalization;
D   = blasso.D; %added

L = size(D,2); %added

A   = blasso.A;

%colex = gencolex(fc); %added

cvx_solver    mosek
cvx_precision high

switch solver
    case 'primal'
        T = zeros(m.^2);
        Psi = zeros(m.^2,L);
        
        cvx_begin sdp %quiet
        variable R(prod(m)+1, prod(m)+1) hermitian
        variable z( n )                  complex
        variable u( prod(2*m-1) )        complex
        variable t(1)                    complex
        dual variable Q
        Q: R >= 0
        % constraints
        for k1=-m(1)+1:m(1)-1
            %id1 = (sum(colex,2) == k1); %added
            
            for k2=-m(2)+1:m(2)-1
                %id2 = (sum(colex,2) == k2); %added
                
                Th = kron( toepl(k2), toepl(k1) );
                T  = T  +  u( k1+m(1)-1 + (2*m(1)-1)*( k2+m(2)-1 ) + 1 ) * Th;
                
                %for l=1:L
                %    Psilk = (D(:,l) .* colex(id))
                %    Psi(:,:,l) = Psi(:,:,l) + u( k1+m(1)-1 + (2*m(1)-1)*( k2+m(2)-1 ) + 1 ) * Psilk; %added
                %end
            end
        end
        R(1:prod(m), 1:prod(m)) == T   ;
        R(I, prod(m)+1)         == z(:);
        R(prod(m)+1, prod(m)+1) == t   ;
        
        minimize( C0 * ( 1/2*real(t) + 1/2*real(R(1,1)) + 1/2/la * gam * norm( y - A(z), 'fro' )^2 ) );
        cvx_end
        
        p = 1/la * (y - A(z));
        
        
    case 'dual' % TODO
        cvx_begin sdp
        variable Q (prod(m)+1,prod(m)+1) hermitian
        variable Q1(prod(m),  prod(m)  ) hermitian
        variable p( n )                  complex
        dual variable R
        R: Q >= 0
        % constraints
        for k1=-m(1)+1:m(1)-1
            for k2=-m(2)+1:m(2)+1
                Th = kron( toepl(k2), toepl(k1) );
                Th(:)' * Q1(:) == double( (k1==0) & (k2==0) );
            end
        end
        
        if lambda==0
            maximize( real(p'*y) )
        else
        end
        cvx_end
        
    otherwise
        error('Unknown solver')
end


end

