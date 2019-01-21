function [F,f0] = fobj(ell,y,lambda,rho,gam,A)
%FOBJ_SR Computes toeplitz-penalized sdp-blasso objective
%   F = FOBJ_SR(fc,y,lambda,rho,gam,A) outputs
%       F  - objective, with low-rank storage: F(U) = f(U*U')
%       f0 - normalization constant, f0 = F(0)
%
%   The objective reads
%   f(R,z,t) = t + tr(R)/m + 1/(2*la)*||y-Az||_H^2 + 1/(2*rho)*||R-Toep(R)||^2
%

m  = 2*ell + 1;
f0 = 2*lambda / gam / norm(y,'fro')^2;


z     = @(U) U(1:end-1,:) * U(end,:)';
normT = @(T) sum( Dnumel1(m) .* abs( T ).^2 );

% TV norm term
F_tv      = @(U) 1/2        * ( norm( U(1:end-1,:), 'fro' )^2 / m + norm( U(end,:), 'fro' )^2 );

% Data-fitting term
F_datafit = @(U) 1/2/lambda * gam * norm( y - A( z(U) ), 'fro' )^2;

% Toeplitz penalization term
F_Tpen1   = @(U) 1/2/rho    * norm( U(1:end-1,:)' * U(1:end-1,:), 'fro' )^2;
F_Tpen2   = @(U) -1/2/rho   * normT( Tproj1( U(1:end-1,:) ) );
F_Tpen    = @(U) F_Tpen1(U) + F_Tpen2(U);



F  = @(U) f0       * ( F_tv(U) + F_datafit(U) + F_Tpen(U) );
%F0 = @(U) C0 * ( F_tv(U) + F_datafit(U) );


%p = @(U) 1/lambda * ( y - A(z) );

end

