function [A,As,S,Amat,Asmat] = approximationOperator(fc,model)
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
%Xf = FreqMesh(fc);
Xf = -fc:fc;

kernel = model.kernel;

switch model.fop
    case 'convolution'
        w = compute_weights(fc,kernel);
        
        A  = @(z) w       .* z;
        As = @(p) conj(w) .* p;
        
        S = @(z) abs(w).^2 .* z;
        
        Amat  = diag(w);
        Asmat = Amat';
        
    case 's-convolution'
        grid  = model.grid;
        
        w  = compute_weights(fc,kernel);
        Xg = compute_grid(fc,grid);
        L  = getoptions(grid,'size',64);
        
        C = sqrt(1/L);
        
        if 2*fc >= L
            error('grid is insufficiently large wrt fc (TODO)');
        end
        
        Amat  = exp( 2i*pi * Xg(:) * Xf(:)' ) * diag(w(:));
        Asmat = C * Amat';
        Smat = Asmat*Amat;
        
        %A  = @(z) Amat * ( w(:) .* z(:) );
        A  = @(z) exp(2i*pi * Xg(:) * Xf(:)') * (w.*z);
        %As = @(p) conj(w) .* (Asmat * p(:));
        As = @(p) (Asmat * p(:));
        
        S = @(z) Smat * z(:);
        
    case 'foveation'
        grid = model.grid;
        
        Xg  = compute_grid(fc,grid);
        L   = length(Xg);
        
        sig = getoptions(kernel, 'cov', @(x)(.001 + .05*abs(x)));
        
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

        switch kernel.mode
            case 1
                phij  = @(x,s) exp( -(s - x).^2 ./ (2*sig(x).^2) );
                phij1 = @(x,s) phij(x-1,s) + phij(x,s) + phij(x+1,s);

                [X,G] = meshgrid(t,Xg);
                B     = phij1(X,G);
                F     = 1/N * fftshift(fft(B,[],2),2);
                Fc    = F( : , (N/2 + 1 - fc) : (N/2 + 1 + fc) );

                Amat  = conj(Fc);
                Asmat = sqrt(1/L) * Amat';

            case 2
                em = sqrt(2*pi) * sig(Xg(:)) .* exp(-2i*pi*Xg(:)*Xf(:)');
                Amat = em .* exp( -2*pi^2*sig(Xg(:)).^2 .* Xf(:)'.^2 );
                Amat = conj(Amat);
                Asmat = sqrt(1/L) * Amat';
        end
        
        A  = @(z) Amat  * z(:);
        As = @(p) Asmat * p(:);
        
        Smat = Asmat * Amat;
        S = @(z) Smat * z(:);
            
    otherwise
        error('TODO');
end


end


function [weights] = compute_weights(fc,kernel)
%COMPUTE_KERNEL_WEIGHTS compute transfer function of filter

d = numel(fc);

switch kernel.type
    case 'Dirichlet'
        n = 2*fc+1;
        weights = ones(n, 1) / n; % normalization l1
        %weights = ones(n,1);
        
    case 'Gaussian'
        % Fourier coefficients of Gaussian are computed in closed-form
        sig = getoptions(kernel, 'cov', .05);
        
        Xf      = FreqMesh(fc);
        weights = ( sqrt(2*pi) * sig )^d * exp( -2*pi^2*sig^2 * Xf(:)'.^2 )';
        
    otherwise
        error('Unknown kernel');
end

end


function [Xg] = compute_grid(fc,grid)
%COMPUTE_SUBSAMPLING_GRID 

switch grid.shape
    case 'lattice'
        L = getoptions(grid, 'size', 64);
        
        q  = ceil( (2*fc+1) / L );
        t  = 1; %TODO
        Lq = L*q;
        %Xg = linspace(0, (Lq-1)*t, Lq) / Lq;
        Xg = linspace(0, (Lq-1)*t, L) / Lq;
        Xg = linspace(0, L-1, L) / L;
        
    otherwise
        error('TODO');
end

end

