function [TX] = Tprod1(T,X)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

d = 1;

s = size(X,1);

pad_size = [ 2*s - 1, size(X,2) ];
pad_X    = Pad1(X, pad_size);

Cprod = ifft( fft( T ) .* fft( pad_X ) );
TX    = Restr1( Cprod, s);


end

