function [g,G0] = fgrad(fc,y,lambda,rho,gam,A,AS)
%FGRAD_SR Computes the gradient of toeplitz-penalized sdp-blasso objective
%   [g,G0] = FGRAD_SR(ell,y,lambda,rho,gam,A,AS) outputs:
%       g  - gradient, with low-rank storage: g(U,h) = < Gradf(U*U'), h >
%       G0 -
%
%   The objective reads
%   f(R,z,t) = t + tr(R)/m + 1/(2*la)*||y-Az||_H^2 + 1/(2*rho)*||R-Toep(R)||^2
%
%   Inputs:
%       fc      - cutoff frequency (R is of size (2*fc+1))
%       y       - measurements
%       lambda  - blasso regularization parameter
%       rho     - toeplitz penalization parameter
%       gam     - constant for the Hilbert space norm
%       A, AS   - approximation operator


n = 2*fc + 1;

ASy = AS(y);
f0  = 2*lambda / gam / norm(y, 'fro')^2; % objective normalization
G0 = 1/(4*prod(n)) + 1/(2*lambda^2) * sum( abs(ASy).^2 ) + 1/4;

z  = @(U) U(1:end-1,:) * U(end,:)';
dz = @(U) 1/2 * ( -AS(y) + AS( A( z(U) ) ) );

% TV norm term
g_tv      = @(   h) 1/2      * [ h(1:end-1,:) / n ; h(end,:) ];

% Data fitting term
g_datafit = @(U, h) 1/lambda * [ dz(U) * h(end,:) ; dz(U)' * h(1:end-1,:) ];

% Toeplitz penalization term
g_Tpen1   = @(U, h) 1/rho    * U(1:end-1,:) * ( U(1:end-1,:)' * h(1:end-1,:) );
g_Tpen2   = @(U, h) -1/rho   * Tprod1( Tproj1( U(1:end-1,:) ), h(1:end-1,:) );
g_Tpen    = @(U, h) [ g_Tpen1(U,h) + g_Tpen2(U,h) ; zeros(1, size(h,2)) ];

%test = @(U,h) zeros(size(h));

g = @(U,h) f0 * ( g_tv(h) + g_datafit(U,h) + g_Tpen(U,h) );
%g = @(U,h) C0 * test(U,h);

%g0 = @(U,h) F0 * ( g_tv(h) + g_datafit(U,h) );

end

