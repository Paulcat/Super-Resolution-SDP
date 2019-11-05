function [N] = Dnumel1(n)
%DNUMEL1 Computes Card{0 <= i,j <= n ; i-j = k}, for k = 0,...,n,-n,...,0

d = 1;

%pad_ones = Pad1( ones( [n(1:d), 1] ), [2*n(1:d)-1, 1] );
P1 = Pad1( ones(n,1) );

N = ifftshift( ifft( fft( P1 ).^2 ) );

end

