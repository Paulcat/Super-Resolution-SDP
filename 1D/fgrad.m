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
%   ELL specifies the size of the variable U, ie U has sizee (2ELL+1) x r.
%   H has size (2ELL+1) x l.
%
%   See also FOBJ, FGRADU

debug = 0;

m = 2*ell + 1;

Asy = As(y);
f0  = 2*la / gam / norm(y, 'fro')^2; % objective normalization
G0 = 1/(4*prod(m)) + 1/(2*la^2) * sum( abs(Asy).^2 ) + 1/4;

% helpers
z  = @(U) U(1:m,:) * U(m+1,:)';
Dz = @(U) 1/2 * (-Asy + AsA( z(U) ));
%Dz = @(U) 1/2 * (z(U) - y);

% gradient
G_bl1 = @(U,h) 1/2/m * h(1:m,:) + 1/la * Dz(U) * h(m+1,:);
G_bl2 = @(U,h) 1/2 * h(m+1,:) + 1/la * Dz(U)' * h(1:m,:);
G_Tp = @(U1,h1) 1/rho * ( U1*(U1'*h1) - Tprod1(Tproj1(U1), h1) );

G = @(U,h) f0 * [G_bl1(U,h) + G_Tp(U(1:m,:),h(1:m,:)); G_bl2(U,h)];
%G = @(U,h) f0 * [1/la*Dz(U)*h(m+1,:); 1/la*Dz(U)'*h(1:m,:)];
%G = @(U,h) f0 * [1/rho * U(1:m,:)*(U(1:m,:)'*h(1:m,:)); zeros(1,size(h,2))];
%G = @(U,h) f0 * [1/rho * Tprod1(Tproj1(U(1:m,:)),h(1:m,:)); zeros(1,size(h,2))];
%G = @(U,h) f0 * [G_Tp(U(1:m,:),h(1:m,:)); zeros(1,size(h,2))];
%G = @(U,h) f0 * [G_bl1(U,h); G_bl2(U,h)];

end


%

% 
% Dz = @(U) 1/2 * ( -Asy + AS( A( z(U) ) ) );
%
% % TV norm term
% g_tv      = @(   h) 1/2      * [ h(1:end-1,:) / m ; h(end,:) ];
% 
% % Data fitting term
% g_datafit = @(U, h) 1/la * [ Dz(U) * h(end,:) ; Dz(U)' * h(1:end-1,:) ];
% 
% % Toeplitz penalization term
% g_Tpen1   = @(U, h) 1/rho    * U(1:end-1,:) * ( U(1:end-1,:)' * h(1:end-1,:) );
% g_Tpen2   = @(U, h) -1/rho   * Tprod1( Tproj1( U(1:end-1,:) ), h(1:end-1,:) );
% g_Tpen    = @(U, h) [ g_Tpen1(U,h) + g_Tpen2(U,h) ; zeros(1, size(h,2)) ];
% 
% %test = @(U,h) zeros(size(h));
% 
% g = @(U,h) f0 * ( g_tv(h) + g_datafit(U,h) + g_Tpen(U,h) );
% %g = @(U,h) C0 * test(U,h);
% 
% %g0 = @(U,h) F0 * ( g_tv(h) + g_datafit(U,h) );
