function coeffs = ls_coeffs(ell,y,la,rho,gam,As,AsA)
%LSCOEFFS_SR Computes line-search coefficients for toeplitz-penalized
%sp-blasso objective
%   [c11,c22,c12,c1,c2] = LSCOEFFS_SR(fc,U,v,blasso) returns the
%   coefficients such that
%
%       F( [a*U, b*v] ) = c11*a^2 + c22*b^2 + c12*a*b + c1*a + c2*b + c0
%
%   which are then used in the line-search step of FFW


Asy = As(y);
C0 = 2*la/gam*norm(y,'fro')^2; % TODO

m = 2*ell+1;


z = @(U) U(1:m,:) * U(m+1,:)';
tau = @(U) norm(U(m+1,:), 'fro')^2;
T = @(U) Tproj1(U(1:m,:));

dotAz = @(U,V) real( z(U)' * AsA(z(V)) );
dotR = @(U,V) sum( abs(U(1:m,:)'*V(1:m,:)).^2, 'all' );
dotT = @(U,V) sum( Dnumel1(m) .* conj(T(U)) .* T(V) );
normT = @(U) sum( Dnumel1(m) .* abs(T(U)).^2 );


% U1   = U(1:end-1,:);
% zU   = U1 * U(end,:)';
% tauU = norm( U(end,:), 'fro')^2;
% v1   = v(1:end-1,:);
% zv   = v1 * v(end,:)';
% tauv = norm( v(end,:), 'fro')^2;


% c11 = C0/2/la * normH( A(zU) )^2     +    C0/2/rho * ...
%     ( norm(U1'*U1, 'fro')^2 - sum( Dnumel1(n) .* abs( Tproj1(U1) ).^2 ) );
c11 = @(U) C0/2/la * dotAz(U,U) + C0/2/rho * ( dotR(U,U) - normT(U) );


% c22 = C0/2/la * normH( A(zv) )^2     +    C0/2/rho * ...
%     ( norm(v1'*v1, 'fro')^2 - sum( Dnumel1(n) .* abs( Tproj1(v1) ).^2 ) );
c22 = @(v) C0/2/la * dotAz(v,v) + C0/2/rho * ( dotR(v,v) - normT(v) );

% c12 = C0/la * real( dotpH( zU, As( A(zv) ) ) ) + C0/rho * ( norm(U1'*v1, 'fro')^2 ...
%     - real( sum( Dnumel1(n) .* conj( Tproj1(U1) ) .* Tproj1(v1) ) ) );
c12 = @(U,v) C0/la*dotAz(U,v) + C0/rho * ( dotR(U,v) - real( dotT(U,v) ) );

%c1 = C0 * ( 1/2*(tauU + norm(U1,'fro')^2/n) - 1/la * real( dotpH( y, A(zU) ) ) );
c1 = @(U) C0/2 * tau(U) + C0/m * norm(U(1:m,:),'fro')^2 - C0/la * real(Asy'*z(U));

%c2 = C0 * ( 1/2*(tauv + norm(v1,'fro')^2/n) - 1/la * real( dotpH( y, A(zv) ) ) );
c2 = @(v) C0/2 * tau(v) + C0/m * norm(v(1:m,:),'fro')^2 - C0/la * real(Asy'*z(v));

coeffs = @(U,v) [c11(U), c22(v), c12(U,v), c1(U), c2(v)];

end

