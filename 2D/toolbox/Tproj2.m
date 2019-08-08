function TU = Tproj2(ell,U)
%TPROJ2 Projection on 2D generalized Toeplitz matrices
%   TU = TPROJ2(U) returns the "diagonal" values of the generalized
%   Toeplitz matrix closest to UU*
%
%   U1 is a matrix of size PROD(M) x r, with M = 2ELL+1, and TU is an array
%   of size (2M(1)-1) x ... x (2M(d)-1) containing all distinct diagonal
%   values, sorted in a fft-friendly order
%
%   Works only if ELL(1) = ... = ELL(d)
%
%   See also TPROD2

d = length(ell); % dimension

% if sum(diff(ell)) ~= 0
%     error('TODO')
% end

r = size(U,2);
m = 2*ell+1;

%U = reshape(U, [m,r]);
%PU   = Pad2(ell, U, [2*m-1,r]);
PU = Pad2(ell, U, 2*m-1);

TU = sum( ifft2( abs( fft2( PU ) ).^2 ) ./ Dnumel2(m), d+1 );


end

