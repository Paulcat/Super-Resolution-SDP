function [TU] = Tproj2(U1)
%TPROJ2 Projection on 2D generalized Toeplitz matrices
%   TU = TPROJ2(U) is an array containing the diagonal values of the
%   generalized-Toeplitz projection of UU*
%
%   The input U must be tensor shaped, of size (n1 x n2 x r). The resulting
%   array TU is of size ( (2*n1-1) x (2*n2-1) ).
%   Actually, the code works only if n1 = n2.
%
%   TODO: ordering??

d = 2; % dimension

n = size(U1);


% if sum(diff(n(1:d))) ~= 0
%     error('TODO')
% end

pad_size = [ 2 * n(1:d) - 1, size(U1,d+1) ];
pad_U1   = Pad2(U1, pad_size);

TU = sum( ifft2( abs( fft2( pad_U1 ) ).^2 ) ./ Dnumel2(n), d+1 );


end

