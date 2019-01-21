function [U2] = TenToMat(U)
%TENTOMAT Reshape tensor as matrix
%   U2 = TENTOMAT(U) -- if U is of size (n1 x n2 x r), TENTOMAT(U) returns 
%   a matrix of size ((n1*n2) x r)

d  = 2;
U2 = reshape( U, numel(U)/size(U,d+1), [] );
%U2 = reshape( U, prod(n), [] );


end

