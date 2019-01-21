function [c11, c22, c12, c1, c2] = lsCoeffs(fc,U,v,blasso)
%LSCOEFFS_SR Computes line-search coefficients for toeplitz-penalized
%sp-blasso objective
%   [c11,c22,c12,c1,c2] = LSCOEFFS_SR(fc,U,v,blasso) returns the
%   coefficients such that
%
%       F( [a*U, b*v] ) = c11*a^2 + c22*b^2 + c12*a*b + c1*a + c2*b + c0
%
%   which are then used in the line-search step of FFW


y   = blasso.y;
la  = blasso.lambda;
rho = blasso.rho;
C0  = blasso.f0;
gam = blasso.ga;

A = blasso.A;
AS = blasso.AS;

dotpH = @(x,y) gam       * (x'*y);
normH = @(x)   sqrt(gam) * norm(x,'fro');

n    = 2*fc+1;
U1   = U(1:end-1,:);
zU   = U1 * U(end,:)';
tauU = norm( U(end,:), 'fro')^2;
v1   = v(1:end-1,:);
zv   = v1 * v(end,:)';
tauv = norm( v(end,:), 'fro')^2;


c11 = C0/2/la * normH( A(zU) )^2     +    C0/2/rho * ...
    ( norm(U1'*U1, 'fro')^2 - sum( Dnumel1(n) .* abs( Tproj1(U1) ).^2 ) );


c22 = C0/2/la * normH( A(zv) )^2     +    C0/2/rho * ...
    ( norm(v1'*v1, 'fro')^2 - sum( Dnumel1(n) .* abs( Tproj1(v1) ).^2 ) );

c12 = C0/la * real( dotpH( zU, AS( A(zv) ) ) ) + C0/rho * ( norm(U1'*v1, 'fro')^2 ...
    - real( sum( Dnumel1(n) .* conj( Tproj1(U1) ) .* Tproj1(v1) ) ) );

c1 = C0 * ( 1/2*(tauU + norm(U1,'fro')^2/n) - 1/la * real( dotpH( y, A(zU) ) ) );

c2 = C0 * ( 1/2*(tauv + norm(v1,'fro')^2/n) - 1/la * real( dotpH( y, A(zv) ) ) );


end
