function [x,mass_atom] = harmouch(mm,nvar,d)
% This is a MATLAB implementation of 
%
% Algorithm 4.1 in
% Harmouch, Khalil, Mourrain, Structured low rank decomposition of multivariate Hankel matrices (2017)
%
% INPUT
% mm: moment matrix (columns indexed by 1 z z^2 z^3 ... and rows indexed
% by 1 conj(z) conj(z^2) conj(z^3) ... in the univariate case)
% nvar: number of variables
% d: order of the moment matrix
%
% OUPUT
% x: atoms
% mass_atom: weight of each atom
%
% Cedric Josz
% LAAS CNRS, Toulouse, August 2017

if rem(d,2) == 0
    d1 = d/2;
    d2 = d/2;
else
    d1 = d/2 - .5;
    d2 = d/2 + .5;
end

% HACK pour gérer le cas d'une matrice nulle
if norm(mm,'fro') < 1e-16
    x = [];
    mass_atom = [];
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate binomial coefficient table
% used to locate a variable from its powers
bin = zeros(nvar,d+2);
bin(1,:) = 0:d+1;
for k = 2:nvar
    bin(k,:) = cumsum(bin(k-1,:));
end

% Generate column vectors of powers, arranged in a matrix
pow = genpow(nvar+1,d);
pow = pow(:,2:nvar+1);

% disp('Moment matrix:')
% mm;

% svd
% disp('Singular value decomposition:')
mmm = mm(1:bin(nvar,d1+2),1:bin(nvar,d2+1));
[U,S,V] = svd(mmm);
mmm - U*S*V';

% rank
eps = 1e-1;
SS = S;
S = diag(S);
r = length(S);
drop = find(S(2:r) ./ S(1:r-1) < eps);
if ~isempty(drop)
    rankM = drop(1);
else
    rankM = r;
end
S = SS;
nb = rankM;

% shift operator
% disp('Multiplication matrices:')
V = conj(V);
T = cell(1,nvar);

for k = 1:nvar
    ext = [];
    for i = 1:bin(nvar,d2+1)
        powj = pow(i,:) + pow(k+1,:);
        j = 1;
        for kk = 1:nvar
            j = j + bin(nvar+1-kk,1+sum(powj(kk:nvar)));
        end
        ext = [ ext j ];
    end
    Ur = U'; Ur = Ur(1:nb,:);
    Vr = conj(V); Vr = Vr(:,1:nb);
    Sr = S(1:nb,1:nb);
    T{k} = Sr^(-1)*Ur*mm(1:bin(nvar,d1+2),ext)*Vr;
    T{k};
end

% disp('Commutativity:')
% T{1}*T{2}-T{2}*T{1}
% T{1}'*T{1}-T{1}*T{1}'
% T{2}'*T{2}-T{2}*T{2}'
% T{1}'*T{2}-T{2}*T{1}'

% Compute common eigenvalues of multiplication matrices

% Random combinations of multiplication matrices
% disp('Choose random coefficients for simultaneous diagonalization:')
if nvar == 1
    coef = 1;
else
    coef = 2*rand(nvar,1)-1;
end
% coef = coef / sum(coef)
% coef = 2*rand(nvar,1)-1+1i*(2*rand(nvar,1)-1)
M = zeros(nb);
for i = 1:nvar
    M = M + coef(i)*T{i};
end;
% disp('Random combination of the shift operators adds up to:')
% M;
%M*M'-M'*M
% disp('Spectral decomposition:')

% ATTENTION RAJOUTS (HACK) DE PAUL
if sum(sum(isinf(M)))==0 && sum(sum(isnan(M)))==0   % HACK!!!
    [Q,D] = eig(M);
else
    Q = zeros(nb);      % HACK!!
end

%[Q,D] = orderschur(M)
%[V,DD] = eig(M)
% if norm(diag(D) - diag(DD)) > 1e-1
%    disp('The ordering of the eigenvectors is wrong.') 
% end
% disp('The shift operators are unitarily similar to:')
for i = 1:nvar
    Q'*T{i}*Q;
end
x = zeros(nvar,1,nb);
mass_atom = zeros(nb,1);
for i = 1:nb
    for j = 1:nvar
        x(j,1,i) = Q(:,i)'*T{j}*Q(:,i);
    vv= [];
    for k = 1:bin(nvar,d2+1)
       vv = [ vv prod((x(:,1,i).').^(pow(k,:))) ];
    end
    e1 = zeros(1,size(mmm,1));
    e1(1,1) = 1;
    mass_atom(i,1) = e1*mmm*Vr*Q(:,i)/(vv*Vr*Q(:,i));
    end
end

end

