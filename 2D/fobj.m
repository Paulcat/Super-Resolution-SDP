function [F,f0] = fobj(ell,y,la,rho,gam,As,AsA)
%FOBJ objective function
%   [F,F0] = FOBJ(ELL,Y,LA,RHO,GAM,AS,ASA) computes the objective F(R,Z,T):
%
%       T + tr(R)/N + 1/(2*LA)*||Y-A(Z)||_H^2 + 1/(2*RHO)*||R-Toep(R)||^2
%
%   so that F(U) = f(UU',Z,T) (T can be removed ?)
%
%   Everything is expressed in terms of A* (AS) and A*A (ASA).
%   GAM is a constant for the Hilbert space norm.
%
%   ELL is a vector specifying the size of the variable U, ie U has size
%   PROD(2ELL+1) x r.
%
%   See also FGRAD

debug = 0;

d  = 2;
m  = 2*ell + 1;
N  = prod(m);
f0 = 2*la / gam / norm(y, 'fro')^2;
Asy = As(y);

% helpers
fro2  = @(x) norm(x,'fro')^2;
normT = @(T) sum( Dnumel2(m) .* abs(T).^2, [1 d] );
z     = @(U) U(1:N,:) * U(N+1,:)';
R     = @(U) U(1:N,:)' * U(1:N,:);

% objective
F_bl = @(U) 1/2/N * fro2(U(1:N,:)) + 1/2 * fro2(U(N+1,:)) + ...
            1/2/la * ( gam*fro2(y) - 2*real(Asy'*z(U)) + real(z(U)'*AsA(z(U))) );
F_Tp = @(U) 1/2/rho * ( fro2(R(U)) - normT( Tproj2(ell,U(1:N,:)) ) );

F = @(U) f0 * ( F_bl(U) + F_Tp(U) );

if debug
    F = @(U) f0 * F_bl(U);
    %F = @(U) f0 * F_Tpen(U);
    %F = @(U) 1/2 * normT( Tproj2(ell,U) );
end

end

%%%%%%%%%%%%%%%%%%%%%

%y = reshape(y,m);

%reshT = @(x) reshape(x, [m, size(x,2)]);
%z = @(U) reshT( U(1:N,:) * U(N+1,:)' );

% F_blasso = @(U) 1/2/N * fro2( U(1:N,:) ) + 1/2 * fro2(U(end,:)) + gam/2/la * fro2( y - A(z(U)) );
