function [TX] = Tprod1(T,X)
%TPROD1 Toeplitz-vector product in 1D
%   TX = TPROD1(T,X) returns the product between the Toeplitz matrix
%   Toep(T) and the matrix X.
%
%   T is an array containing the diagonal values of the Toeplitz matrix

global nfft
nfft = nfft + 3;

%[m,r] = size(X);

%PX       = Pad1(X, [2*m-1, r]);
PX = Pad1(X);

Cprod = ifft( fft( T ) .* fft( PX ) );
TX    = Restr1(Cprod);


end

