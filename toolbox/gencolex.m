function [mom_colex] = gencolex(fc)
%GENCOLEX(fc) Returns moment orders (powers) sorted following the
%colexicographic order
%
%   fc: vector (size d) of maximum orders (in absolute value), in each
%   dimension.
%
%   NOTE: analogous of genpow, which does the same for the graded
%   lexicographic order.

d = numel(fc);

switch d
    case 1
        mom_colex = (-fc:fc)';
    
    case 2
        [k1,k2] = meshgrid(-fc(1):fc(1) ,-fc(2):fc(2));
        k1 = k1';
        k2 = k2';
        mom_colex = [k1(:),k2(:)];
        
    case 3
        [k1,k2,k3] = meshgrid(-fc(1):fc(1), -fc(2):fc(2), -fc(3):fc(3));
        k1 = permute(k1,[2 1 3]);
        k2 = permute(k2,[2 1 3]);
        mom_colex = [k1(:),k2(:),k3(:)];
        
    otherwise
        error('TODO');
end

end

