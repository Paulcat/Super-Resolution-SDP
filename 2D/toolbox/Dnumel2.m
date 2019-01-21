function [N] = Dnumel2(n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

d = 2; % dimension

pad_ones = Pad2( ones( [n(1:d), 1] ), [2*n(1:d)-1, 1] );

N = ifftshift( ifft2( fft2( pad_ones ).^2 ) );

end

