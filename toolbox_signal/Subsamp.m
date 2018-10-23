function [SX] = Subsamp(X,p)
%SUBSAMP Downsample matrix X of factor p = p1 x ... x pd.
%   X has size n1 x ... x nd x r.

% probably, taking into account the last dimension r is not useful, as it
% will likely always be 1

dim = numel(p);

switch dim
    case 1
        SX = X(1:p:end,:);
    case 2
        SX = X(1:p(1):end,1:p(2):end,:);
    case 3
        SX = X(1:p(1):end,1:p(2):end,1:p(3):end,:);
    otherwise
        error('TODO')
end

end

