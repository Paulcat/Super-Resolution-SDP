function TU = Tproj1(U)
%TPROJ1 Toeplitz projection (in 1D)
%   TU = TPROJ1(U) is an array containing the diagonal values of the
%   Toeplitz projection of U*U'.
%
%   If U is of size (N x R), TU is of size ((2*N-1) x R).
%   Diagonal values are sorted as 0, 1, ..., n-1, -n+1, ..., -1.

global nfft;
nfft = nfft + 2;

%d = 1;

%[m,r] = size(U);
n = size(U,1);

%PU       = Pad1(U, [2*m-1,r]);
PU = Pad1(U);

%TU = sum( ifft( abs( fft( PU ) ).^2 ) ./ Dnumel1(m), d+1 );
TU = sum( ifft( abs(fft(PU)).^2 ) ./ Dnumel1(n), 2 );

end

