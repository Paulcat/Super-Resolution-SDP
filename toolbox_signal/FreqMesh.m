function [Xf,Yf,Zf] = FreqMesh(fc)
%FREQMESH Summary of this function goes here
%   Detailed explanation goes here

dim = numel(fc);

switch dim
    case 1
        Xf = (-fc:fc)';
    case 2
        [Xf,Yf] = meshgrid(-fc(1):fc(1), -fc(2):fc(2));
    case 3
        [Xf,Yf,Zf] = meshgrid(-fc(1):fc(1), -fc(2):fc(2), -fc(3):fc(3));
    otherwise
end


end

