function [UX] = Upsamp2(n,X,p)
%UPSAMP2 2D upsampling operator
%   UX = UPSAMP2(N,X,P), with X of size (PROD(N) x r), returns a "tensor"
%   of size (N(1)*P(1) x N(2)*P(2) x r) corresponding to the
%   zero-upsampling of X by a factor P.
%
%   See also SUBSAMP2. Called by APPROXIMATIONOPERATOR.

%d = 2;

r = size(X,2);

X = reshape(X,[n,r]);
UX = zeros( [p.*n,r] );
UX(1:p(1):end, 1:p(2):end, :) = X;

% s  = size(X);
% s2 = s(1:d);
% 
% UX = zeros( [ p.*s2, size(X,d+1) ] );
% UX(1:p(1):end, 1:p(2):end, :) = X;


end

