function [collapsing] = check_collapsing(U,ell)
%CHECK_COLLAPSING simple check for hierarchy collapsing
%   CHECK_COLLAPSING(U,l) checks the flatness condition from [Curto and
%   Fialkow, 1996] on matrix U*U'. U should be of size (prod(l)+1).

N1 = prod( 2*ell + 1 ); % main block size
N2 = prod( 2*(ell-1) + 1 ); % sub-block size

U1 = U(1:N1,:);
U2 = U(1:N2,:);

sv1 = svd(U1).^2; % (non-zero) singular values of U1*U1' (other svalues are necessarily zeros)
sv2 = svd(U2).^2; % same thing for U2*U2'

nsv1 = sv1 ./ max(sv1);
nsv2 = sv2 ./ max(sv2);

r1 = sum( nsv1 > 1e-3 );
r2 = sum( nsv2 > 1e-3 );

collapsing = (r1==r2);


end

