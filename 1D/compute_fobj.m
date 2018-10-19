function [F,F0] = compute_fobj(fc,y,lambda,rho,gam,A)
%COMPUTE_FOBJ compute approximate blasso objective
%   F = COMPUTE_FOBJ(fc,y,lambda,rho,gam,A) returns F

n = 2*fc + 1;

C0 = 2*lambda / gam / norm(y,'fro')^2;


z     = @(U) U(1:end-1,:) * U(end,:)';
normT = @(T) sum( Dnumel1(n) .* abs( T ).^2 );

F_tv      = @(U) 1/2        * ( norm( U(1:end-1,:), 'fro' )^2 / n + norm( U(end,:), 'fro' )^2 );
F_datafit = @(U) 1/2/lambda * gam * norm( y - A( z(U) ), 'fro' )^2;
F_Tpen1   = @(U) 1/2/rho    * norm( U(1:end-1,:)' * U(1:end-1,:), 'fro' )^2;
F_Tpen2   = @(U) -1/2/rho   * normT( Tproj1( U(1:end-1,:) ) );
F_Tpen    = @(U) F_Tpen1(U) + F_Tpen2(U);

F  = @(U) C0       * ( F_tv(U) + F_datafit(U) + F_Tpen(U) );
%F = @(U) C0 * 0;

F0 = @(U) C0 * ( F_tv(U) + F_datafit(U) );


%p = @(U) 1/lambda * ( y - A(z) );

end

