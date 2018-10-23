function [PX] = Pad(X,padSize)
%PAD Padding operator
%   If X is of size n1 x ... x nd x r, Pad(X) returns zero-padded tensor of 
%   size p1 x ... x pd x r as specified by padsize. d must be strictly
%   greater than 1, r must be explicitly specified (even if it is 1);

PX = zeros(padSize);

switch numel(padSize)
    case 1
        error('Dimension cannot be 0');
    case 2
        PX(1:size(X,1),:) = X;
    case 3
        PX(1:size(X,1),1:size(X,2),:) = X;
    case 4
        PX(1:size(X,1),1:size(X,2),1:size(X,3),:) = X;
    otherwise
        error('TODO')
end


end

