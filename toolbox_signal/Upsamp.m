function [UX] = Upsamp(X,p)
%UPSAMP Upsample matrix X of factor p = p1 x ... x pd, by inserting zeros.
%   X has size n1 x ... x nd x r.

% probably, taking into account the last dimension r is not useful, as it
% will likely always be 1

dim = numel(p);

s  = size(X);
sd = s(1:dim);

UX = zeros([p.*sd, size(X,dim+1)]);

switch dim
    case 1
        UX(1:p:end,:) = X;
    case 2
        UX(1:p(1):end,1:p(2):end,:) = X;
    case 3
        UX(1:p(1):end,1:p(2):end,1:p(3):end,:) = X;
    otherwise
        error('TODO');
end


end

