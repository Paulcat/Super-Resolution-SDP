function [R,z] = perform_approx_blasso_sdp(fc,blasso,options)
%PERFORM_APPROX_BLASSO_SDP solve the toeplitz-penalized sdp-blasso
%   Solve the following program:
%
%   min t + tr(R)/m + 1/(2*la)*||y-Az||^2 + 1/(2*rho)*||R - Toep(R)||^2
%   s.t.    |R  |z|
%           |___|_| >= 0
%           |z* |t|
%
%   Inputs: 
%       fc    - cutoff frequency (R is of size (2*fc+1))
%       blasso - blasso.y, etc... specifies the quantity composing the
%       objective

n = 2*fc + 1;
m = n;


y   = blasso.y;
la  = blasso.lambda;
rho = blasso.rho;
gam = blasso.ga;
C0  = blasso.f0;
A   = blasso.A;

cvx_solver SeDuMi
cvx_precision high


cvx_begin sdp quiet
variable R(m+1,m+1) hermitian
variable R1(m,m)    hermitian toeplitz
variable z(n)       complex
dual variable Q
Q : R >= 0
%
R(1:n,m+1) == z
if la==0
    error('TODO');
else
    minimize( C0 * ( 1/2*(R1(1,1) + R(end,end)) + ...
        1/2/la  * gam * norm(y - A(z), 'fro')^2 + ...
        1/2/rho * norm(R(1:m,1:m) - R1,'fro')^2 ) );
end

cvx_end


end

