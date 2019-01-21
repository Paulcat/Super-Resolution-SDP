function [A,AS,Amat,ASmat] = approximationOperator(fc,model)
%APPROXIMATIONOPERATOR Compute the Fourier approximation matrix A
%   [A,AS] = APPROXIMATIONOPERATOR(fc,model) returns the approximation
%   operator A (and its adjoint AS) of the forward operator specified in
%   model, in the Fourier basis (e_k), -fc <= k <= fc, where e_k: x ->
%   exp(-2i*pi*k*x)
%
%   fc:     cutoff frequency
%   model:  model.fop       - forward operator type (convolution |
%           subsampled-convolution | foveation)
%           model.kernel    - kernel type (Gaussian | Dirichlet)
%           model.kparam    - specifies the variance for Gaussian kernel
%           model.gsize     - size of the subsampling grid , for subsampled
%           measurements (only for subsampled-convolution or foveation)
%
%   For more info on the approximation matrix, see
%   Catala, P., Duval, V. and Peyré, G., A low-rank approach to
%   off-the-grid sparse super-resolution

% n = 2*fc + 1;
Xf = FreqMesh(fc);

switch model.fop
    case 'convolution'
        w = compute_kernel_weights(fc,model);
        
        A  = @(z) w       .* z;
        AS = @(p) conj(w) .* p;
        
        Amat  = diag(w);
        ASmat = Amat';
        
    case 'subsampled-convolution'
        w  = compute_kernel_weights  (fc,model);
        Xg = compute_subsampling_grid(fc,model);
        L  = length(Xg);
        
        Amat  = exp( 2i*pi * Xg(:) * Xf(:)' ) * diag(w(:));
        ASmat = sqrt(1/L) * Amat';
        
        A  = @(z) Amat * ( w(:) .* z(:) );
        AS = @(p) conj(w(:)) .* (ASmat * p(:));
        
    case 'foveation'
        Xg  = compute_subsampling_grid(fc,model);
        L   = length(Xg);
        
        sig = getoptions(model, 'kparam', @(x)(0.001 + .05*abs(x)));
        
        N = 2048;
        t = (0:(N-1))'/N;

% for debug purpose ---
%         Amat = zeros(L,2*fc+1);
        
%         phij  = @(X,sj) exp( -(sj - X(:)).^2 ./ (2*sig(X(:)).^2) );
%         phij1 = @(X,sj) phij(X-1,sj) + phij(X,sj) + phij(X+1,sj);
%         
%         for j=1:L
%             uj = phij1(t,Xg(j));
%             %plot(X,uj), axis tight, drawnow;
%             Fj = 1/N * fftshift(fft(uj)); % fft approximation (in the spirit of Riemann sum)
%             
%             %figure(2), plot(X,uj), drawnow;
%             %figure(1), plot(X,real(Fj)), drawnow;
%             Fj = Fj( (N/2 + 1 - fc) : (N/2 + 1 + fc) );
%             % normalize?
%             
%             Amat(j,:) = conj(Fj);
%         end
% ----------------------
        
        phij  = @(X,Xg) exp( -(Xg - X).^2 ./ (2*sig(X).^2) );
        phij1 = @(X,Xg) phij(X-1,Xg) + phij(X,Xg) + phij(X+1,Xg);
        
        [X,G] = meshgrid(t,Xg);
        B     = phij1(X,G);
        F     = 1/N * fftshift(fft(B,[],2),2);
        Fc    = F( : , (N/2 + 1 - fc) : (N/2 + 1 + fc) );
        
        Amat  = conj(Fc);
        ASmat = sqrt(1/L) * Amat';
        
        A  = @(z) Amat  * z(:);
        AS = @(p) ASmat * p(:);
            
    otherwise
        error('TODO');
end


end


function [weights] = compute_kernel_weights(fc,model)
%COMPUTE_KERNEL_WEIGHTS compute transfer function of filter

d = numel(fc);

switch model.kernel
    case 'Dirichlet'
        n = 2*fc+1;
        weights = ones(n, 1) / n; % normalization l1
        
    case 'Gaussian'
        % Fourier coefficients of Gaussian are computed in closed-form
        sig = getoptions(model, 'kparam', .05);
        
        Xf      = FreqMesh(fc);
        weights = ( sqrt(2*pi) * sig )^d * exp( -2*pi^2*sig^2 * Xf(:)'.^2 )';
        
    otherwise
        error('Unknown kernel');
end

end


function [Xg] = compute_subsampling_grid(fc,model)
%COMPUTE_SUBSAMPLING_GRID 

switch model.grid
    case 'lattice'
        L = getoptions(model, 'gsize', 64);
        
        q  = ceil( (2*fc+1) / L );
        t  = 1; %TODO
        Lq = L*q;
        %Xg = linspace(0, (Lq-1)*t, Lq) / Lq;
        Xg = linspace(0, (Lq-1)*t, L) / Lq;
        
    otherwise
        error('TODO');
end

end

