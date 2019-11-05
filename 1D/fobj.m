function [F,f0] = fobj(n,y,la,rho,gam,As,AsA)
%FOBJ objective function (in 1D)
%   [F,F0] = FOBJ(N,Y,LA,RHO,GAM,AS,ASA) computes the objective f(R,Z,T):
%
%       T + tr(R)/M + 1/(2*LA)*||Y-A(Z)||_H^2 + 1/(2*RHO)*||R-Toep(R)||^2
%
%   so that F(U) = f(UU',Z,T)
%
%   The objective is expressed in terms of A* (AS) and A*A (ASA).
%   GAM is a constant for the Hilbert space norm.
%
%   N specifies the main size of U
%
%   See also FGRAD

debug = 0;

%m  = 2*ell + 1;

% helpers
fro2  = @(x) norm(x,'fro')^2;
normT = @(T) sum( Dnumel1(n) .* abs(T).^2 );
z     = @(U) U(1:n,:) * U(n+1,:)';
R     = @(U) U(1:n,:)' * U(1:n,:);


f0  = 2*la / gam / fro2(y);
Asy = As(y);


% objective
F_bl = @(U) 1/2/n * fro2(U(1:n,:)) + 1/2 * fro2(U(n+1,:)) + ...
    1/2/la * ( gam*fro2(y) - 2*real(Asy'*z(U)) + real(z(U)'*AsA(z(U))) );

F_Tp = @(U) 1/2/rho * ( fro2(R(U)) - normT( Tproj1(U(1:n,:)) ) );

F = @(U) f0 * (F_bl(U) + F_Tp(U));
%F = @(U) f0 * 1/2/rho * fro2(R(U));
%F = @(U) f0 * 1/2/rho * normT(Tproj1(U(1:m,:)));
%F = @(U) f0 * F_Tp(U);
%F = @(U) f0 * F_bl(U);

end


%

% % TV norm term
% F_tv      = @(U) 1/2        * ( norm( U(1:end-1,:), 'fro' )^2 / m + norm( U(end,:), 'fro' )^2 );
% 
% % Data-fitting term
% F_datafit = @(U) 1/2/la * gam * norm( y - A( z(U) ), 'fro' )^2;
% 
% % Toeplitz penalization term
% F_Tpen1   = @(U) 1/2/rho    * norm( U(1:end-1,:)' * U(1:end-1,:), 'fro' )^2;
% F_Tpen2   = @(U) -1/2/rho   * normT( Tproj1( U(1:end-1,:) ) );
% F_Tpen    = @(U) F_Tpen1(U) + F_Tpen2(U);
%
% F  = @(U) f0       * ( F_tv(U) + F_datafit(U) + F_Tpen(U) );
% %F0 = @(U) C0 * ( F_tv(U) + F_datafit(U) );
%
% %p = @(U) 1/lambda * ( y - A(z) );
