function [Gpre] = fgradU(U,ell,y,la,rho,gam,As,AsA)
%FGRADU_SR Compute gradient of toep-penalized sdp-blasso objective at U*U'
%   g_pre = FGRADU_SR(U,fc,y,lambda,rho,gam,A,AS) is the functional such
%   that
%       g_pre(h) = < Gradf(U*U'), h >
%   In particular, the Toeplitz projection of U*U' is pre-computed

m = 2*ell + 1;

f0  = 2*la / gam / norm(y,'fro')^2; % objective normalization
Asy = As(y);

% pre-computations
U1  = U(1:m,:);
z   = U1 * U(m+1,:)';
Dz  = 1/2 * ( -Asy + AsA(z) );
T   = Tproj1(U1);

% pre-computed gradient operator
G_bl1 = @(h) 1/2/m * h(1:m,:) + 1/la * Dz * h(m+1,:);
G_bl2 = @(h) 1/2 * h(m+1,:) + 1/la * Dz' * h(1:m,:);
G_Tp  = @(h1) 1/rho * ( U1*(U1'*h1) - Tprod1(T,h1) );

Gpre = @(h) f0 * [G_bl1(h) + G_Tp(h(1:m,:)); G_bl2(h)];
%Gpre = @(h) f0 * [G_bl1(h); zeros(1,size(h,2))];

end

%
% g1 = @(h) [ 1/2 * h(1:end-1,:) / m        + 1/rho * U1 * (U1' * h(1:end-1,:)); ...
%             1/2 * h(end,:) ];
%         
% g2 = @(h) [ 1/lambda * Dz  * h(end,    :) - 1/rho * Tprod1(T, h(1:end-1,:)); ...
%             1/lambda * Dz' * h(1:end-1,:) ];
% 
%         
% g_pre  = @(h) f0 * ( g1 (h) + g2 (h) );
