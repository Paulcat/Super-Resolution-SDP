function [U] = MatToTen(U2,n)
%MATTOTEN Reshape matrix as tensor (2D)
%   U = MATTOTEN(U2,n)
%
%   Inputs:
%       n = [n1,n2]
%       U2 of size (prod(n) x r)
%
%   Outputs a tensor-shaped matrix of size (n1 x n2 x r)

U = reshape(U2, [n, size(U2,2)]);


end

