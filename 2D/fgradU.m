function [Gpre] = fgradU(U,ell,y,la,rho,gam,As,AsA)
%FGRADU gradient of the objective with precomputed parts
%   GPRE = FGRADU(U,ZETA,ELL,FC,Y,LA,RHO,GAM,A,AS) computes the gradient at
%   U of the objective f(R,Z,T)
%   
%       T + tr(R)/N + 1/(2*LA)*||Y-A(Z)||_H^2 + 1/(2*RHO)*||R-Toep(R)||^2
%
%   so that GPRE(H) = gradf(UU',Z,T) * H. In particular, the Toeplitz
%   projection of UU' is pre-computed.
%
%   GAM is a constant for the Hilbert space norm.
%
%   ELL is a vector specifying the size of the variable H, ie H has size
%   PROD(2ELL+1) x l.
%
%   See also FGRAD

d = 2;
m = 2*ell + 1;
N = prod (m);

f0  = 2*la / gam / norm(y, 'fro')^2; % objective normalization
Asy = As(y);

% helpers
flat = @(x) x(:);
%reshT = @(x) reshape(x, [m,size(x,2)]);
%reshM = @(t) reshape(t, [N,size(t,d+1)]);

% pre-computations
U1  = U(1:N,:);
z   = U(1:N,:)*U(N+1,:)';
T   = Tproj2(ell,U1);
Dz  = 1/2 * flat( -Asy + AsA(z) );

% "pre-computed" gradient operator
G_bl1 = @(h) 1/2/N * h(1:N,:) + 1/la * Dz * h(end,:);
G_bl2 = @(h) 1/2 * h(end,:) + 1/la * Dz' * h(1:N,:);
G_Tpen = @(h1) 1/rho * ( U1*(U1'*h1) - Tprod2(ell,T,h1) );

Gpre = @(h) f0 * [G_bl1(h) + G_Tpen(h(1:N,:)); G_bl2(h)];

end

