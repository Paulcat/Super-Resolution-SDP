function [M1] = Kron_to_Lexico(M)
%KRON_TO_CANONICAL Rearrange the moment matrix from [Dumitrescu] in the
%canonical basis (as in [Lasserre])
%   Useful only for the multivariate case (in the univariate case, the two
%   orderings are identical)
%   This code is designed only for the 2D case!
%   M should be ordered following (1 z2 z2^2 ...) \otimes (1 z1 z1^2 ...)
%   M1 is ordered as (1 z1 z2 z1^2 z1z2 z2^2 .... z2^kmax)

d = 2;  % dimension
n = sqrt(size(M,1));    % n = 2f_c + 1
kmax = n-1;  % kmax = 2f_c, maximal moment order

TM = reshape(M, [n n n n]); % Tensor form of the "Kronecker" moment matrix

s = nchoosek(d+kmax,d); % size of the canonical basis
pow = genpow(d+1,kmax);
pow = pow(:,2:d+1) + 1; % +1 because indices begin at 1 in matlab

M1 = zeros(s);
for i=1:s
    for j=1:s
        i1 = pow(i,1); i2 = pow(i,2);
        j1 = pow(j,1); j2 = pow(j,2);
        M1(i,j) = TM(i1,i2,j1,j2);
    end
end
% NB: not all elements in TM are taken: many are discarded... (diamond
% shape, ..?)

end

