function [eVec,eVal,infos] = perform_LMO_step(options,Gprod,v0)
%PERFORM_LMO_STEP linear minimization oracle step in main algorithm fwsdp
%   [eVec,eVal,infos] = PERFORM_LMO_STEP(options,gprod,v0) computes the 
%   minimal eigenvalue and a corresponding eigenvector of the gradient, 
%   with Power Iterations.
%
%   gprod(h) returns < Grad, h >


maxIter = getoptions(options,'maxIter',1000);

% Initialization
v = v0;
cos_th = 0.5; % cosin of angle between current and previous vectors
niter1 = 0;
niter2 = 0;

%while abs( abs(cos_th) - 1 ) > 1e-24   &&   niter1 < maxIter
while niter1 < maxIter
    x = Gprod(v);
    
    % normalize
    x = x ./ norm(x,'fro');
    
    % update
    cos_th = x'*v;
    v      = x;
    niter1 = niter1 + 1;
end

eVal1 = v'*Gprod(v); % largest (in absolute value) eigen value of the gradient

if eVal1 < 0 % then eVal1 is the lowest eigenvalue
    eVal = eVal1;
    eVec = v;

    
    

else % then eVal1 is the largest eigenvalue, and the lowest can be find using PI on (G - eVal1*Id)
    v = v0;
    cos_th = 0.5;
    
    %while abs( abs(cos_th) - 1 ) > 1e-20    &&   niter2 < maxIter
    while niter2 < maxIter
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

