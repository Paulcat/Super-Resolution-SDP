function [X] = Restr1(PX,n)
%RESTR1 Truncation operator in 1D
%   X = Restr1(PX,n), returns the submatrix X of PX of size n = n1 x r. r
%   must be explicitely specified, even when r=1.


X = PX(1:n,:);


end

