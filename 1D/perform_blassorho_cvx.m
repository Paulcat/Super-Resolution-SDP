function [R,z] = perform_blassorho_cvx(fc,blasso,options)
%PERFORM_BLASSORHO_CVX solve the approximate sdp blasso
%   Solve the SDP program associated to the Toeplitz-penalized blasso

n = 2*fc + 1;
m = n;


y   = blasso.y;
la  = blasso.lambda;
rho = blasso.rho;
gam = blasso.ga;
C0  = blasso.normalization;
A   = blasso.A;

cvx_solver mosek
cvx_precision high


cvx_begin sdp
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

