function [g,G0] = fgrad(ell,y,lambda,rho,gam,A,AS)
%FGRAD_SR Computes the gradient of toeplitz-penalized sdp-blasso objective
%   [g,G0] = COMPUTE_FGRAD(fc,y,lambda,rh,gam,A,AS) outputs:
%       g  - gradient, with low-rank storage: g(U,h) = < Gradf(U*U'), h >
%       G0 -
%
%   The objective reads
%   f(R,z,t) = t + tr(R)/m + 1/(2*la)*||y-Az||_H^2 + 1/(2*rho)*||R-Toep(R)||^2
%
%   Inputs:
%       ell    - order in the Lasserre hierarchy (variable R is of size
%       (2*ell+1)^2. ell should be greater than the cutoff frequency fc
%       y      - measurements
%       lambda - blasso regularization parameter
%       rho    - toeplitz penalization parameter
%       gam    - constant for the Hilbert space norm
%       A, AS  - approximation operator (and adjoint)

d   = 2;
m   = 2*ell + 1;
N   = prod(m);

f0  = 2*lambda / gam / norm(y, 'fro')^2; % obkective normalization
ASy = AS(y);
G0  = 1/(4*N) + 1/(2*lambda^2) * sum( abs(ASy).^2 ) + 1/4;

flat = @(x) x(:);

% Tensor-shaped variables z and U
z    = @(U) MatToTen( U(1:end-1,:) * U(end,:)', m );
dz   = @(U) 1/2 *  (  -ASy   +  AS( A( z(U) ) )  );
Tens = @(U) MatToTen( U(1:end-1,:), m );


% TV norm term
g_tv      = @(   h) 1/2      * [ h(1:end-1,:) / N ; h(end,:) ];

% Data fitting term
g_datafit = @(U, h) 1/lambda * [ flat(dz(U)) * h(end,:) ; flat(dz(U))' * h(1:end-1,:) ];

% Toeplitz penalization term
g_Tpen1   = @(U, h) 1/rho    * U(1:end-1,:) * ( U(1:end-1,:)' * h(1:end-1,:) );
g_Tpen2   = @(U, h) -1/rho   * TenToMat( Tprod2( Tproj2( Tens(U) ), Tens(h) ) );
g_Tpen    = @(U, h) [ g_Tpen1(U,h) + g_Tpen2(U,h); zeros(1, size(h,2)) ];


g = @(U,h) f0 * ( g_tv(h) + g_datafit(U,h) + g_Tpen(U,h) );
%g = @(U,h) g_datafit(U,h);

%g = @(U) 2 * C0 * ( g_tv(U) + g_datafit(U,U) + g_Tpen(U,U) );
%G = @(U) 2 * g(U,U);

end

