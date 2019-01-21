function [SX] = Subsamp2(X,p)
%SUBSAMP2 Dowsample X of factor p
%   SUBSAMP2(X,p), with X of size n1 x n2 x r, and p = p1 x p2
%
%   See also UPSAMP2


SX = X(1:p(1):end, 1:p(2):end, :);

end

