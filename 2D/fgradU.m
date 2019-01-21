function [g_pre] = fgradU(U,fc,y,lambda,rho,gam,A,AS)
%FGRADU_SR Compute gradient of toep-penalized sdp-blasso objective at U*U'
%   g_pre = FGRADU_SR(U,fc,y,lambda,rho,gam,A,AS) is the functional such
%   that
%       g_pre(h) = < Gradf(U*U'), h >
%   In particular, the Toeplitz projection of U*U' is pre-computed

d = 2;
n = 2*fc + 1;
N = prod (n);


Tens  = @(U) MatToTen( U(1:end-1,:), n );

f0  = 2*lambda / gam / norm(y, 'fro')^2; % objective normalization
ASy = AS(y);
U1  = U(1:end-1,:);
z   = MatToTen( U1 * U(end,:)', n );
dz  = 1/2 * ( -ASy  +  AS( A(z) ) );
T   = Tproj2( Tens(U) );


g1 = @(h) [ 1/2 * h(1:end-1,:)/N             + 1/rho * U1 * (U1' * h(1:end-1,:)); ...
            1/2 * h(end,:) ];
g2 = @(h) [ 1/lambda * dz(:)  * h(end,    :) - 1/rho * TenToMat( Tprod2( T, Tens(h) ) ) ; ...
            1/lambda * dz(:)' * h(1:end-1,:) ];

g_pre = @(h) f0 * ( g1(h) + g2(h) );

end

