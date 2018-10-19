function [aa,bb] = perform_linesearch_step(fc,Ut,vt,blasso)
%PERFORM_LINESEARCH_STEP linesearch step in main algorithm fwsdp
%   [mu,nu] = PERFORM_LINESEARCH_STEP(fc,U,v,blasso) returns the solution
%   of min_(mu,nu) f( mu*(U*U') + nu*(v*v') )
%
%   See also FWSDP


[cxx, cyy, cxy, cx, cy] = compute_ls_coefficients(fc,Ut,vt,blasso);

denom = 4*cxx*cyy - cxy^2;
num1  = cxy*cy - 2*cyy*cx;
num2  = cxy*cx - 2*cxx*cy;

if denom == 0
    if cxx == 0
        mu = 0;
        nu = min(1, max(0, -cy/(2*cyy)));
    else    % HACK?
        mu = min(1, max(0, -cx/(2*cxx)));
        nu = 0;
    end
else
    mu = num1/denom;
    nu = num2/denom;
    
    if ~(mu>=0) || ~(nu>=0) || ~(mu+nu<=1)
        % optimal point on the segment a+b=1?
        mu1 = max(0, min(1, -(cxy - 2*cyy + cx - cy)/(2*(cxx - cxy + cyy))));
        nu1 = 1 - mu1;
        f1 = cxx*mu1^2 + cxy*mu1*nu1 + cyy*nu1^2 + cx*mu1 + cy*nu1;
        
        % optimal point on the ordinate segment [0,1]?
        mu2 = 0;
        nu2 = max(0, min(1, -cy/(2*cyy)));
        f2 = cyy*nu2^2 + cy*nu2;
        
        % optimal point on the absciss segment [0,1]?
        mu3 = max(0, min(1, -cx/(2*cxx)));
        nu3 = 0;
        f3 = cxx*mu3^2 + cx*mu3;
        
        [~,I] = min([f1,f2,f3]);
        ind = zeros(3,1); ind(I) = 1;
        mu = [mu1 mu2 mu3] * ind;
        nu = [nu1 nu2 nu3] * ind;
    end
end




if denom == 0
    if cxx == 0
        a = 0;
        b = min(1, max(0, -cy/(2*cyy)));    %HACK!
    else                %HACK!
        a = min(1, max(0, -cx/(2*cxx)));
        b = 0;
    end
else
    a = num1/denom;
    b = num2/denom;

    if a+b>=1
        a = max(0, min(1, -(cxy - 2*cyy + cx - cy)/(2*(cxx - cxy + cyy))));
        b = 1-a;
    else
        if a<=0 && b<=0
            a = 0;
            b = 0;
        elseif a<=0
            a = 0;
            b = max(0, min(1, -cy/(2*cyy)));
        elseif b<=0
            a = max(0, min(1, -cx/(2*cxx)));
            b = 0;
        end
    end
end


% CVX check
A = [cxx, cxy/2; cxy/2, cyy];
B = [cx, cy];
S = sqrtm(A);
%cvx_solver mosek
cvx_precision high
cvx_begin quiet
variable X(2) nonnegative
% 'simplex' constraints
sum(X) <= 1
%
minimize( norm(S*X,'fro')^2 + B*X )
cvx_end

aa = X(1);
bb = X(2);

% abs(aa-mu)/abs(aa)
% abs(bb-nu)/abs(nu)

if abs(a-mu) > 1e-12 || abs(b-nu) > 1e-12
   warning('pb precision ls?');
end



end