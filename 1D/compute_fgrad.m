function [g,g0] = compute_fgrad(fc,y,lambda,rho,gam,A,AS)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


n = 2*fc + 1;

C0  = 2*lambda / gam / norm(y, 'fro')^2;
ASy = AS(y);

z  = @(U) U(1:end-1,:) * U(end,:)';
dz = @(U) 1/2 * ( -AS(y) + AS( A( z(U) ) ) );



g_tv      = @(   h) 1/2      * [ h(1:end-1,:) / n ; h(end,:) ];
g_datafit = @(U, h) 1/lambda * [ dz(U) * h(end,:) ; dz(U)' * h(1:end-1,:) ];
g_Tpen1   = @(U, h) 1/rho    * U(1:end-1,:) * ( U(1:end-1,:)' * h(1:end-1,:) );
g_Tpen2   = @(U, h) -1/rho   * Tprod1( Tproj1( U(1:end-1,:) ), h(1:end-1,:) );
g_Tpen    = @(U, h) [ g_Tpen1(U,h) + g_Tpen2(U,h) ; zeros(1, size(h,2)) ];

%test = @(U,h) zeros(size(h));

g = @(U,h) C0 * ( g_tv(h) + g_datafit(U,h) + g_Tpen(U,h) );
%g = @(U,h) C0 * test(U,h);

g0 = @(U,h) C0 * ( g_tv(h) + g_datafit(U,h) );

end

