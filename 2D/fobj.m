function [F,f0] = fobj(fc,y,lambda,rho,gam,A)
%FOBJ_SR Computes toeplitz-penalized sdp-blasso objective
%   F = FOBJ_SR(fc,y,lambda,rho,gam,A) outputs
%       F  - objective, with low-rank storage: F(U) = f(U*U')
%       f0 - normalization constant, f0 = F(0)
%
%   The objective reads
%   f(R,z,t) = t + tr(R)/m + 1/(2*la)*||y-Az||_H^2 + 1/(2*rho)*||R-Toep(R)||^2
%

d  = 2;
n  = 2*fc + 1;
N  = prod(n);
f0 = 2*lambda / gam / norm(y, 'fro')^2;

flat  = @(x) x(:);

% Tensor-shaped variables z and U
z     = @(U) MatToTen( U(1:end-1,:) * U(end,:)', n );
Tens  = @(U) MatToTen( U(1:end-1,:), n );
normT = @(T) sum( flat( Dnumel2(n) .* abs( T ).^2 ) );


F_tv      = @(U) 1/2        * ( norm( U(1:end-1,:), 'fro' )^2 / N + norm( U(end,:), 'fro' )^2 );
F_datafit = @(U) 1/2/lambda * gam * norm( y - A(z(U)), 'fro' )^2;
F_Tpen    = @(U) 1/2/rho    * ( norm( U(1:end-1,:)' * U(1:end-1,:), 'fro' )^2 - normT( Tproj2( Tens(U) ) ) );


F = @(U) f0       * ( F_tv(U) + F_datafit(U) + F_Tpen(U) );
%F = @(U) F_tv(U);
%p = @(U) 1/lambda * ( y - A(z(U)) );

end

