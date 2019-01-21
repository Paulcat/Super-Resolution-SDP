function [TX] = Tprod1(T,X)
%TPROD1 Toeplitz-vector product in 1D
%   TX = TPROD1(T,X) returns the product between the Toeplitz matrix
%   Toep(T) and the matrix X.
%
%   T is an array containing the diagonal values of the Toeplitz matrix

global nfft
nfft = nfft + 3;

d = 1;

s = size(X,1);

pad_size = [ 2*s - 1, size(X,2) ];
PX       = Pad1(X, pad_size);

Cprod = ifft( fft( T ) .* fft( PX ) );
TX    = Restr1( Cprod, s);


end

