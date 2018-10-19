function bool = isToeplitz(fc,U1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

tol = 1e-6;

n = 2*fc+1;

T     = Tproj1(U1);
normT = sum( Dnumel1(n) .* abs(T).^2 );
normU = norm( U1'*U1, 'fro')^2;

bool = (normU - normT) < tol;


end

