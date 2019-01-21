function [SX] = Subsamp1(X,s)
%SUBSAMP1 Subsampling operator in 1D
%   SX = SUBSAMP1(X,s) subsamples of X by a factor s

SX = X(1:s:end,:);

end

