function [TU] = Tproj1(U)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

d = 1;

n = size(U);

pad_size = [ 2 * n(1:d) - 1, size(U,d+1) ];
pad_U    = Pad(U, pad_size);

TU = sum( ifft( abs( fft( pad_U ) ).^2 ) ./ Dnumel1(n), d+1 );

end

