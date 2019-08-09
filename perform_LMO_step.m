function [eVec,eVal,infos] = perform_LMO_step(options,Gprod,v0)
%PERFORM_LMO_STEP linear minimization oracle step in ffw algorithm
%   [eVec,eVal,infos] = PERFORM_LMO_STEP(options,gprod,v0) computes the 
%   minimal eigenvalue (eVal) and a corresponding eigenvector (eVec) of the 
%   gradient, with Power Iterations.
%
%   Inputs:
%       Gprod   - returns < Grad(U), h >, at current U
%       v0      - initial vector for power iterarions
%       options - options.maxIter
%                 options.lmoTol


maxIter = getoptions(options,'maxIter',1000);
tol     = getoptions(options,'lmoTol',1e-10);

% Initialization
v = v0;
cos_th = 0.5; % cosin of angle between current and previous vectors
niter1 = 0;
niter2 = 0;

while abs( abs(cos_th) - 1 ) > tol   &&   niter1 < maxIter
%while niter1 < maxIter
    x = Gprod(v);
    
    % normalize
    x = x ./ norm(x,'fro');
    
    % update
    cos_th = x'*v;
    v      = x;
    niter1 = niter1 + 1;
end

eVal1 = v'*Gprod(v); % largest (in absolute value) eigen value of the gradient

if eVal1 < 0 % eVal1 is the lowest eigenvalue
    eVal = eVal1;
    eVec = v;

    
    

else % eVal1 is the largest eigenvalue: the lowest is found using PI on (G - eVal1*Id)
    v = v0;
    cos_th = 0.5;
    
    while abs( abs(cos_th) - 1 ) > tol    &&   niter2 < maxIter
    %while niter2 < maxIter
        x = Gprod(v) - eVal1*v;
        
        % normalize
        x = x ./ norm(x,'fro');
        
        % update
        cos_th = x'*v;
        v      = x;
        niter2 = niter2 + 1;
    end
    
    eVal = v'*Gprod(v);
    eVec = v;

    
    
end




infos.niter = niter1 + niter2;
infos.precision = abs( abs(cos_th) - 1 );

% Eigen space error
err = norm(Gprod(v) - eVal*v,'fro') / abs(eVal);
infos.error = err;

end

