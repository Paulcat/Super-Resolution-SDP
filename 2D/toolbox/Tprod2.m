function [TX] = Tprod2(T,X)
%TPROD2 Toeplitz-vector product in 2D
%   TX = TPROD2(T,X) returns the product between the generalized Toeplitz
%   matrix Toep(T) and the matrix X.
%
%   T and X must be tensor-shaped, ie X is of size n1 x n2 x r, and T is of
%   size (2*n1-1) x (2*n2-1).

d = 2; % dimension

s = size(X);

pad_size = [ 2 * s(1:d) - 1, size(X,d+1) ];
pad_X    = Pad2(X, pad_size);


Cprod = ifft2( fft2(T) .* fft2( pad_X ) );
TX    = Restr2( Cprod, s );


end

