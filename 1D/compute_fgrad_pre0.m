function [g_pre0] = compute_fgrad_pre0(U,fc,y,lambda,gam,A,AS)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

n = 2*fc + 1;

C0  = 2*lambda / gam / norm(y,'fro')^2;
ASy = AS(y);
U1  = U(1:end-1,:);
z   = U1 * U(end,:)';
dz  = 1/2 * ( -ASy + AS( A( z ) ) );

g1 = @(h) [ 1/2 * h(1:end-1,:) / n ; 1/2 * h(end,:) ];
g2 = @(h) [ 1/lambda * dz  * h(end,    :) ; 1/lambda * dz' * h(1:end-1,:) ];

        
g_pre0  = @(h) C0 * ( g1 (h) + g2 (h) );

end