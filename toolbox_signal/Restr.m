function [X] = Restr(PX,truncSize)
%RESTR Truncation operator
%   If PX is of size p1 x ... x pd x r, Restr(X) returns zero-padded tensor 
%   of size n1 x ... x nd x r as specified by truncsize

switch numel(truncSize)
    case 1
        error('Dimension cannot be 0');
    case 2
        X = PX(1:truncSize,:);
    case 3
        X = PX(1:truncSize(1),1:truncSize(2),:);
    case 4
        X = PX(1:truncSize(1),1:truncSize(2),1:truncSize(3),:);
    otherwise
        error('TODO');
end

end

