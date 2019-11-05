function [PX] = Pad1(X)
%PAD1 Padding to double size (in 1D)
%   PX = PAD1(X) when X is a N x R matrix returns PX = [X;0] such that PX
%   has size (2N-1) x R
%
%   Internal use only
%   see also RESTR1; TPROJ1, TPROD1, DNUMEL1

[n,r] = size(X);
p = 2*n-1;

PX = zeros(p,r);
PX(1:n, :) = X;


end

