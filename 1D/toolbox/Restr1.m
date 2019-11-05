function [X] = Restr1(PX)
%RESTR1 Truncating to half size (in 1D)
%   X = Restr1(PX) when PX is a matrix of size (2N-1) x R returns the upper
%   submatrix of PX of size N x R.
%
%   Internal use only
%   See also PAD1; TPROD1

m = size(PX,1);
n = (m+1)/2;

X = PX(1:n,:);


end

