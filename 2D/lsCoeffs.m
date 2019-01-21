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
ga  = blasso.ga;

A   = blasso.A;
AS  = blasso.AS;

dotpH = @(x,y) ga       * (x'*y);
normH = @(x)   sqrt(ga) * norm(x,'fro');
sum2  = @(x)   sum( sum(x) );
flat  = @(x)   x(:);

n     = 2*fc + 1;
N     = prod(n);
U1    = U(1:end-1,:);
zU    = MatToTen( U1 * U(end,:)', n );
TensU = MatToTen(U1, n);
tauU  = norm( U(end,:), 'fro')^2;
v1    = v(1:end-1,:);
zv    = MatToTen( v1 * v(end,:), n );
Tensv = MatToTen(v1, n);
tauv  = norm( v(end,:), 'fro')^2;

% ** c11 = 2||z(UU')||_H^2 + 1/(2*rho) * ||UU' - Toep(UU')||^2 **
c11 = C0/2/la * normH( A(zU) )^2   +   C0/2/rho * ...
    ( sum2( abs(U1'*U1).^2 ) - sum2( Dnumel2(n).*abs( Tproj2( TensU ) ).^2 ) );

% ** c22 = 2||z(vv')||_H^2 + 1/(2*rho) * ||vv' - Toep(vv')||^2 **
c22 = C0/2/la * normH( A(zv) )^2   +   C0/2/rho * ...
    ( sum2( abs(v1'*v1).^2 ) - sum2( Dnumel2(n).*abs( Tproj2( Tensv ) ).^2 ) );

% ** c12 = 4real<z(U),z(v)>_H + 1/rho * real<UU'-toep(UU'),vv'-toep(vv')> **
        % simplify: <UU'-toep(UU'),vv'-toep(vv')> = <UU',vv'> - <toep(UU'),toep(vv');
c12 = C0/la*real( dotpH( zU(:), flat( AS( A(zv) ) ) ) ) + C0/rho * ( norm(U1'*v1,'fro')^2 ...
    - real( sum2( Dnumel2(n) .* conj( Tproj2( TensU ) ) .* Tproj2( Tensv ) ) ) );

% ** c1 = tau(U) + 1/n*trace(UU'(1:end-1,1:end-1)) + 2/la*real<y,z(U)>_H **
c1 = C0 * (1/2*(tauU + norm(U1,'fro')^2/N) - 1/la*real( dotpH( y(:), flat( A( zU ) ) ) ));

% ** c2 = tau(v) + 1/n*trace(vv'(1:end-1,1:end-1)) + 2/la*real<y,z(v)>_H **
c2 = C0 * (1/2*(tauv + norm(v1,'fro')^2/N) - 1/la*real( dotpH( y(:), flat( A( zv ) ) ) ));



end

