function TX = Tprod2(ell,T,X)
%TPROD2 Toeplitz-vector product in 2D
%   TX = TPROD2(ELL,T,X) returns the product between the generalized
%   Toeplitz matrix Toep(T) and the matrix X.
%
%   X is a matrix/vector of size PROD(M) x r, with M = 2ELL+1, and T is an
%   array of size 2M(1)-1 x ... x 2M(d)-1 containing the diagonal values of
%   the Toeplitz matrix. TX is a matrix/vector of size PROD(M) x r.
%
%   See also TPROJ2

%d = 2; % dimension

r = size(X,2);
m = 2*ell+1;

%X = reshape(X, [m,r]);
%PX = Pad2(ell, X, [2*m-1,r]);
PX = Pad2(ell, X, 2*m-1);

Cprod = ifft2( fft2(T) .* fft2( PX ) );

%TX    = Restr2( Cprod, m );
TX = Cprod(1:m(1),1:m(2),:);
TX = reshape(TX, [prod(m), r]);


end

