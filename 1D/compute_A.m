function [A,AS,w,Amat,ASmat] = compute_A(fc,model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = 2*fc + 1;
Xf = FreqMesh(fc);


switch model.problem
    case 'convolution'
        w = compute_kernel_weights(fc,model);
        
        A  = @(z) w       .* z;
        AS = @(p) conj(w) .* p;
        
        Amat  = eye(n);
        ASmat = eye(n);
        
    case 'subsampled-convolution'
        w  = compute_kernel_weights  (fc,model);
        Xg = compute_subsampling_grid(fc,model);
        L  = length(Xg);
        
        Amat  = exp( 2i*pi * Xg(:) * Xf(:)' ); % here, A is simply F* over the sampling grid
        ASmat = sqrt(1/L) * Amat';
        
        A  = @(z) Amat * ( w(:) .* z(:) );
        AS = @(p) conj(w(:)) .* (ASmat * p(:));
        
    case 'foveation'
%         Xg  = compute_subsampling_grid(fc,model);
%         L   = length(Xg);
        
        Xg = linspace(0,1,512);
        L = length(Xg);
        
        sig = getoptions(model, 'kparam', @(x)(0.001 + abs(x)));
        
        N = 2048;
        %X = [0:N/2, -N/2:-1]' / N; ?
        X = (0:(N-1))'/N;
        Amat = zeros(L,2*fc+1);
        
        for j=1:L
            Xj = X;
            phij  = exp( -(Xj(:)-Xg(j)).^2 ./ (2*sig(Xj(:)).^2) );
            Fphij = 1/N * fftshift(fft(phij)); % fft approximation (in the spirit of Riemann sum)
            Fphij = Fphij( (N/2 + 1 - fc) : (N/2 + 1 + fc) );
            % normalize?
            
            Amat(j,:) = conj(Fphij);
        end
        
        ASmat = sqrt(1/L) * Amat';
        
        A  = @(z) Amat  * z(:);
        AS = @(p) ASmat * p(:);
        
        w = 0; %TODO!        
    otherwise
        error('TODO');
end


end


function [weights] = compute_kernel_weights(fc,model)
%COMPUTE_KERNEL_WEIGHTS compute transfer function of convolution filter
%   In the case of *convolution*, one may describe the observations in the
%   Fourier domain

d = numel(fc);

switch model.kernel
    case 'Dirichlet'
        n = 2*fc+1;
        weights = ones(n, 1) / n; % normalization l1
        
    case 'Gaussian'
        % close form for Fourier coefficients
        sig = getoptions(model, 'kparam', .05);
        
        % specific to 1D
        Xf      = FreqMesh(fc);
        weights = ( sqrt(2*pi) * sig )^d * exp( -2*pi^2*sig^2 * Xf(:)'.^2 )';
        
    otherwise
        error('Unknown kernel');
end

end


function [Xg] = compute_subsampling_grid(fc,model)

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

