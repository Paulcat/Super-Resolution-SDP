function [UX] = Upsamp2(X,p)
%UPSAMP2 Upsample X of factor p (by inserting zeros)
%   UPSAMP2(X,p), with X of size n1 x n2 x r, and p = p1 x p2
%
%   See also SUBSAMP2

d = 2;

s  = size(X);
s2 = s(1:d);

UX = zeros( [ p.*s2, size(X,d+1) ] );
UX(1:p(1):end, 1:p(2):end, :) = X;


end

