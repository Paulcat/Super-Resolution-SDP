function [PX] = Pad2(X,PadSize)
%PAD2 Padding operator (2D)
%   PAD2(X,PADSIZE) (outer) pads X up to size PADSIZE
%
%   X must be passed in tensor form (n1 x n2 x r), and PADSIZE must follow
%   this format (p1 x p2 x r), with r explicit (even if it is 1)
%
%   See also RESTR2

PX = zeros(PadSize);
PX(1:size(X,1), 1:size(X,2), :) = X;


end

