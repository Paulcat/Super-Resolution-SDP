function PX = Pad2(ell,X,p)
%PAD2 Padding operator (2D)
%   PX = PAD2(ELL,X,P), where X is a matrix of size PROD(2*ELL+1) x r,
%   returns a "tensor" of size P(1) x P(2) x r, where each component of
%   size P(1) x P(2) is th result of outer-zero-padding the columns of X,
%   seen as squared matrices.
%
%   I leave the result under its tensor form as PAD2 is excusively used as
%   an intermediary step: see e.g. Tproj2, Tprod2
%
%   See also RESTR2. Called by TPROD2, TPROJ2, DNUMEL2

if sum(p < ell)
    error('error: size after padding is smaller than before')
end

r = size(X,2);
X = reshape(X,[2*ell+1,r]);

PX = zeros([p,r]);
PX(1:2*ell(1)+1, 1:2*ell(2)+1, :) = X;
%PX = reshape(PX,[prod(p),r]);


end

