function [d,f,X,A] = flat_norm(x,D2,y,b,tau)

% flat_norm - compute flat norm
%
%   [d,f,X,A] = flat_norm(x,a,y,b,tau);
%
%   d is the flat norm between 
%       sum_i a_i delta_{x_i} 
%       sum_j b_j delta_{y_j}
%
%   x should be of size (n,d) where d is the dimension.
%
%   _Warning:_ the method is slower in dimension d>1. 
%
%   Copyright (c) 2017 Gabriel Peyre

if nargin<5
    % flatness parameter
    tau = 1;
end
    
d = size(x,2);

X = [x; y];
A = [D2(:);-b(:)];
P = size(X,1);    
    
if d==1
    %%% special 1D case %%%
    [X,I] = sort(X); A = A(I);        
    cvx_solver sdpt3 % SeDuMi %
    cvx_begin sdp quiet
    cvx_precision high;
    variable f(P,1);
    abs(f(2:end)-f(1:end-1)) <= X(2:end)-X(1:end-1);
    norm(f,Inf)<=tau;
    maximize( A'*f )
    cvx_end
    d = A'*f;
    return;
end

%%% Arbitrary dimension case %%%

% pairwise distance matrix.
D2 = sum(X.^2,2);
D = repmat(D2,1,P) + repmat(D2',P,1) - 2*X*X';
D = sqrt(max(D,0));
% all pair of constraints
[J,I] = meshgrid(1:P,1:P);
s = find(I<J); I = I(s); J = J(s);

cvx_solver sdpt3 % SeDuMi %
cvx_begin sdp quiet
cvx_precision high;
variable f(P,1);
abs(f(I)-f(J)) <= D(I+P*(J-1));
norm(f,Inf)<=tau;
maximize( A'*f );
cvx_end
d = A'*f;


end


