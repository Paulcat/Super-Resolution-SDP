function X = Restr2(N,PX,T)
%RESTR Truncation operator (2D)
%   X = RESTR2(N,PX,T), where PX is a matrix of size (PROD(N) x r),  
%
%(outer) truncates X down to size TRUNCSIZE
%
%   PX must be passed in tensor form (p1 x p2 x r), and TRUNCSIZE must be
%   of the form (n1 x n2)
%
%   See also PAD2. Called by TPROD2.
% Useless???

X = PX(1:T(1),1:T(2),:);

end

