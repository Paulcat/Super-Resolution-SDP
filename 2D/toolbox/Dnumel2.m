function D = Dnumel2(n)
%DNUMEL2 Counting diagonal elements
%   D = DNUMEL2(N) counts the number of "diagonal" elements in a
%   2D generalized Toeplitz matrix of size prod(N) x prod(N)

ell = (n-1)/2;

%P1 = Pad2( ones( [n, 1] ), [2*n-1, 1] );
P1 = Pad2(ell, ones(prod(n),1), 2*n-1);

D = ifftshift( ifft2( fft2( P1 ).^2 ) );

end

