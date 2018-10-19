function [N] = Dnumel1(n)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

d = 1;

pad_ones = Pad( ones( [n(1:d), 1] ), [2*n(1:d)-1, 1] );

N = ifftshift( ifft( fft( pad_ones ).^2 ) );

end

