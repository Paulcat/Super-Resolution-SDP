function [X] = Restr2(PX,TruncSize)
%RESTR Truncation operator (2D)
%   RESTR2(X,TRUNCSIZE) (outer) truncates X down to size TRUNCSIZE
%
%   PX must be passed in tensor form (p1 x p2 x r), and TRUNCSIZE must be
%   of the form (n1 x n2)
%
%   See also PAD2


X = PX(1:TruncSize(1),1:TruncSize(2),:);

end

