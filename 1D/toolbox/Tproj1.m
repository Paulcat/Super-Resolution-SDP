function [TU] = Tproj1(U)
%TPROJ1 Projection on Toeplitz matrices
%   TU = TPROJ1(U) is an array containing the diagonal values of the
%   Toeplitz projection of U*U'.
%
%   If U is of size (n x r), TU is of size ((2*n-1) x r). Diagonal values
%   are sorted as 0, -1, ..., -n+1, n-1, ..., 1.

global nfft;
nfft = nfft + 2;

d = 1;

[m,r] = size(U);

PU       = Pad1(U, [2*m-1,r]);

TU = sum( ifft( abs( fft( PU ) ).^2 ) ./ Dnumel1(m), d+1 );

end

