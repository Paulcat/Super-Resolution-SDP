function [G,G0] = fgrad(ell,y,la,rho,gam,As,AsA)
%FGRAD gradient of the objective
%   [G,G0] = FGRAD(ELL,Y,LA,RHO,GAM,AS,ASA) computes the gradient of the
%   objective f(R,Z,T):
%   
%       T + tr(R)/N + 1/(2*LA)*||Y-A(Z)||_H^2 + 1/(2*RHO)*||R-Toep(R)||^2
%
%   so that G(U,H) = gradf(UU',Z,T) * H.
%
%   Everything is expressed in terms of A* (AS) and A*A (ASA). GAM is a 
%   constant for the Hilbert space norm.
%
%   ELL is a vector specifying the size of the variable U, ie U has size
%   (PROD(2ELL+1) + 1) x r. H has size (PROD(2ELL(1)+1) + 1) x l.
%
%   See also FGRADU, FOBJ

debug = 0;

d   = 2;
m   = 2*ell + 1;
N   = prod(m);

f0  = 2*la / gam / norm(y, 'fro')^2; % objective normalization
Asy = As(y);

% helpers
flat = @(x) x(:);
z = @(U) U(1:N,:) * U(N+1,:)';
Dz = @(U) 1/2 * flat(-Asy + AsA( z(U) ) );

% gradient
G_bl1 = @(U,h) 1/2/N * h(1:N,:) + 1/la * Dz(U) * h(end,:);
G_bl2 = @(U,h) 1/2 * h(end,:) + 1/la * Dz(U)' * h(1:N,:);
G_Tp = @(U1,h1) 1/rho * ( U1*(U1'*h1) - Tprod2(ell,Tproj2(ell,U1), h1) );

G = @(U,h) f0 * [G_bl1(U,h) + G_Tp(U(1:N,:),h(1:N,:)); G_bl2(U,h)];

if debug
    G = @(U,h) f0 * [G_bl1(U,h); G_bl2(U,h)];
    %G = @(U,h) f0 * [G_Tpen(U(1:N,:),h(1:N,:)); zeros(1,size(h,2))];
    %G = @(U,h) Tprod2(ell,Tproj2(ell,U),h);
end

% gradient in 0
G0  = 1/(4*N) + 1/(2*la^2) * sum( abs(Asy).^2 ) + 1/4;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = reshape(y,m); 
%reshT = @(x) reshape(x, [m,size(x,2)]);
%z = @(U) reshT( U(1:N,:) * U(N+1,:)' );
%Dz = @(U) 1/2 * flat(-ASy + AS( A( z(U) ) ));