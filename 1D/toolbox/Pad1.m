function [PX] = Pad1(X,p)
%PAD1 Padding operator in 1D
%   PX = PAD1(X,p) returns PX = [X;0] such that PX has size p
%   If X is of size n1 x r, p should be of the form p1 x r (r must be
%   explicitely specified, even if r=1.

PX = zeros(p);
PX(1:size(X,1), :) = X;


end

