function coeffs = ls_coeffs(ell,y,la,rho,gam,As,AsA)
%LS_COEFFS line-search coefficients
%   COEFFS = LS_COEFFS(ELL,Y,LA,RHO,GAM,AS,ASA) returns the coefficients
%   such that
%
%       F( [a*U, b*v] ) = c11*a^2 + c22*b^2 + c12*a*b + c1*a + c2*b + c0
%
%   where F is the Blasso-Toeplitz objective.
%
%   See also FOBJ, PERFORM_LINESEARCH_STEP


Asy = As(y);
C0  = 2*la/gam*norm(y,'fro')^2; % TODO

m     = 2*ell + 1;
N     = prod(m);

z = @(U) U(1:N,:) * U(N+1,:)';
tau  = @(U) norm( U(N+1,:), 'fro')^2;
T = @(U) Tproj2(ell,U(1:N,:));

dotAz = @(U,V) real( z(U)' * AsA(z(V)) );
dotR  = @(U,V) sum( abs(U(1:N,:)'*V(1:N,:)).^2, 'all');
dotT  = @(U,V) sum( Dnumel2(m) .* conj(T(U)) .* T(V), 'all' );
normT = @(U)   sum( Dnumel2(m) .* abs(T(U)).^2, 'all');

% ** c11 = 2||z(UU')||_H^2 + 1/(2*rho) * ||UU' - Toep(UU')||^2 **
c11 = @(U) C0/2/la * dotAz(U,U) + C0/2/rho * ( dotR(U,U) - normT(U) );

% ** c22 = 2||z(vv')||_H^2 + 1/(2*rho) * ||vv' - Toep(vv')||^2 **
c22 = @(v) C0/2/la * dotAz(v,v) + C0/2/rho * ( dotR(v,v) - normT(v) );

% ** c12 = 4real<z(U),z(v)>_H + 1/rho * real<UU'-toep(UU'),vv'-toep(vv')> **
        % simplify: <UU'-toep(UU'),vv'-toep(vv')> = <UU',vv'> - <toep(UU'),toep(vv');
c12 = @(U,v) C0/la*dotAz(U,v) + C0/rho * ( dotR(U,v) - real( dotT(U,v) ) );

% ** c1 = tau(U) + 1/n*trace(UU'(1:end-1,1:end-1)) + 2/la*real<y,z(U)>_H **
c1 = @(U) C0/2*tau(U) + C0/N*norm(U(1:N,:),'fro')^2 - C0/la*real(Asy'*z(U));

% ** c2 = tau(v) + 1/n*trace(vv'(1:end-1,1:end-1)) + 2/la*real<y,z(v)>_H **
c2 = @(v) C0/2*tau(v) + C0/N*norm(v(1:N,:),'fro')^2 - C0/la*real(Asy'*z(v));


coeffs = @(U,v) [c11(U), c22(v), c12(U,v), c1(U), c2(v)];


end

