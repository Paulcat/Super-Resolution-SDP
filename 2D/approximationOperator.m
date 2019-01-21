function [A,AS,Amat,ASmat] = approximationOperator(fc,model)
%APPROXIMATIONOPERATOR Compute the Fourier approximation matrix A
%   [A,AS] = APPROXIMATIONOPERATOR(fc,model) returns the approximation
%   operator A (and its adjoint AS) of the forward operator specified in
%   model, in the Fourier basis (e_k), -fc <= (k) <= fc, where e_k: x ->
%   exp(-2i*pi*<k,x>)
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


n = 2*fc+1;

%order   = fft_to_colex(fc);
%fftsort = @(x) x(order);


[Yf,Xf] = FreqMesh(fc);
%[Xg,Yg] = compute_subsampling_grid(fc,model);
%figure, scatter(Xg(:),Yg(:),'.');


switch model.fop
    case 'convolution'
        w = compute_kernel_weights(fc,model);
        
        A  = @(z) w.*z;
        AS = @(p) conj(w).*p;
        
        Amat  = diag(w(:));
        ASmat = diag(w(:))';
        
    case 'subsampled-convolution'
        w       = compute_kernel_weights  (fc,model);
        [Xg,Yg] = compute_subsampling_grid(fc,model);
        L       = getoptions(model, 'gsize', [64 64]);
        
        Amat  = exp( 2i*pi * (Xg(:)*Xf(:)' + Yg(:)*Yf(:)') ) * diag(w(:));
        ASmat = sqrt(1/prod(L)) * Amat';
        
        if strcmp(model.grid, 'lattice') % fast implementation available
            q  = ceil( (2*fc+1) ./ L );
            t  = 1; % TODO!
            Lq = L .* q;
            eL = exp( -2i*pi* ( fc(1)*Xg + fc(2)*Yg ) );
            
            A  = @(z) Subsamp2( prod(Lq) * (eL .* ifft2( Pad2( w .* z, [Lq, 1] ) ) ), t*q );
            AS = @(p) sqrt(1/prod(L)) * conj(w) .* Restr2( fft2( conj(eL) .* Upsamp2( p, t*q ) ), [n, 1] );
            % prod(Lq) is here to counterbalance ifft2
            
        else % no fast implementation, store the entire matrix
            A  = @(z) reshape( Amat  * ( w(:) .* z(:) )    , [L,1] );
            AS = @(p) reshape( conj(w(:)) .* (ASmat * p(:)), [n,1] );
            
        end
        
    case 'foveation'
        [Xg,Yg] = compute_subsampling_grid(fc,model);
        L       = getoptions(model, 'gsize' , [64 64]);
        sig     = getoptions(model, 'kparam', @(x,y)(0.001 + x.^2 + y.^2) );
        
        
        N = 150; % size of the grid over which we compute the fft
        [Y,X] = meshgrid( (0:(N-1))' / N );
        Amat = zeros(prod(L), prod(2*fc+1));
        
        wb = waitbar(0,'Computing matrix A');
        figure(1), clf;
        
        % Compute covariance matrix
        % decompose multiplication with inverse covariance matrix
        % covariance matrix S is defined as
        %   S = T' * D * T, where
        %
        %   T = |cos(t) -sin(t)| and D = |d1 0|
        %       |sin(t)  cos(t)|         |0 d2|
        % t, d1 and d2 may depend on x
        t    = @(x,y) atan2(y-0.5, x-0.5);
        r    = @(x,y) sqrt( (x-0.5).^2 + (y-0.5).^2 );
        d1 = @(x,y) 1 ./ (.001 + .0015 * r(x,y));
        d2 = @(x,y) 1 ./ (0.015)^2;
        
        Si11 = @(x,y) d1(x,y) .* cos(t(x,y)).^2 + d2(x,y) .* sin(t(x,y)).^2;
        Si22 = @(x,y) d1(x,y) .* sin(t(x,y)).^2 + d2(x,y) .* cos(t(x,y)).^2;
        Si12 = @(x,y) (-d1(x,y) + d2(x,y)) .* cos(t(x,y)) .* sin(t(x,y));
        
        vSiv = @(x,y,sx,sy) Si11(x,y) .* (x-sx).^2 + ...
                            2 * Si12(x,y) .* (x-sx) .* (y-sy) + ...
                            Si22(x,y) .* (y-sy).^2;
                        
        % periodized kernel
        phij1 = @(sx,sy) exp( -1/2 * vSiv(X-1, Y-1, sx, sy) ) + ...
                         exp( -1/2 * vSiv(X-1, Y  , sx, sy) ) + ...
                         exp( -1/2 * vSiv(X-1, Y+1, sx, sy) ) + ...
                         exp( -1/2 * vSiv(X  , Y-1, sx, sy) ) + ...
                         exp( -1/2 * vSiv(X  , Y  , sx, sy) ) + ...
                         exp( -1/2 * vSiv(X  , Y+1, sx, sy) ) + ...
                         exp( -1/2 * vSiv(X+1, Y-1, sx, sy) ) + ...
                         exp( -1/2 * vSiv(X+1, Y  , sx, sy) ) + ...
                         exp( -1/2 * vSiv(X+1, Y+1, sx, sy) );
        
        for j=1:prod(L)
            uj = phij1( Xg(j), Yg(j) );
            %figure(1), surf(X,Y,uj,'LineStyle','none'); view(2); colorbar;
            Fj = 1/N^2 * fftshift( fft2(uj) ); % TODO: in this case, close form
            Fj = Fj( (N/2+1-fc(1)) : (N/2+1+fc(1)), (N/2+1-fc(2)) : (N/2+1+fc(2)) );
            % normalize?
            
            Amat(j,:) = conj(Fj(:));
            
            waitbar(j/prod(L),wb);
        end
        close(wb);
        
        ASmat = sqrt(1/prod(L)) * Amat';
        
        A  = @(z) reshape( Amat  * z(:), [L,1] );
        AS = @(p) reshape( ASmat * p(:), [n,1] );

%     TODO --------------------
%     case 'meg-eeg'
%         [Xg,Yg] = compute_subsampling_grid(fc,model);
%         L       = getoptions(model, 'gsize',  [64 64]);
%         sig     = getoptions(model, 'kparam', 0.02);
%         r       = getoptions(model, 'gparam', .2);
%         %sig = 1;
%         
%         N = 256; % size of the grid over which we compute the fft
%         %[Y,X] = meshgrid( (0:N-1)' / N );
%         [X,Y] = ugrid_disc([.5,.5],r-1/(N+1),[N,N]);
%         Amat = zeros(prod(L), prod(2*fc+1));
%         
%         % normalized singular kernel
%         %I = @(X,Y) 1 ./ sqrt(2*pi*sig^4 * (1 + X.^2 + Y.^2) ./ (1 - X.^2 - Y.^2).^3); %* ((X.^2 + Y.^2) <= 1);
%         I = @(X,Y) 1 ./ sqrt(2*pi*sig^4 * (r + (X-0.5).^2 + (Y-0.5).^2) ./ (r - (X-0.5).^2 - (Y-0.5).^2).^3); %* ((X.^2 + Y.^2) <= 1);
%         
%         %phij = @(X,Y,sjx,sjy) I(X,Y) .* ( 1 ./ (1 + 1/sig^2 * ((sjx - X).^2 + (sjy - Y).^2)) );
%         phij = @(X,Y,sjx,sjy) I(X,Y) .* ( sig^2 ./ ( (sjx - X).^2 + (sjy - Y).^2) );
%         phij1 = @(X,Y,sjx,sjy) phij(X,   Y,   sjx, sjy) + ...
%                                phij(X-1, Y,   sjx, sjy) + ...
%                                phij(X+1, Y,   sjx, sjy) + ...
%                                phij(X,   Y-1, sjx, sjy) + ...
%                                phij(X,   Y+1, sjx, sjy) + ...
%                                phij(X-1, Y-1, sjx, sjy) + ...
%                                phij(X+1, Y+1, sjx, sjy) + ...
%                                phij(X-1, Y+1, sjx, sjy) + ...
%                                phij(X+1, Y-1, sjx, sjy);
%         
%         %
%         
%         %figure(1), clf;
%         %figure(2), clf;
%         for j=1:prod(L)
%             uj = phij1(X,Y,Xg(j),Yg(j));
%             %figure(1), surf(X,Y,uj,'LineStyle','none'); zlim([0 1200]);
%             Fj = 1/N^2 * fftshift( fft2(uj) );
%             %figure(2), surf(X,Y,abs(Fj),'LineStyle','none'); zlim([0 8]); view(2); drawnow;
%             Fj = Fj( (N/2+1-fc(1)) : (N/2+1+fc(1)), (N/2+1-fc(2)) : (N/2+1+fc(2)) );
%             % normalize?
%             
%             Amat(j,:) = conj(Fj(:));
%         end
%         
%         ASmat = sqrt(1/prod(L)) * Amat';
%         
%         A  = @(z) reshape( Amat  * z(:), [L,1] );
%         AS = @(p) reshape( ASmat * p(:), [n,1] );
%     ---------------------------


    otherwise
        error('TODO');
          
end


end


function [weights] = compute_kernel_weights(fc,model)

d = numel(fc);

switch model.kernel
    case 'Dirichlet'
        n = 2*fc+1;
        %weights = ones(prod(n), 1) / prod(n); % normalization l1
        weights = ones(n) / prod(n); % normalization l1
        
    case 'Gaussian'
        sig = getoptions(model, 'kparam', .05);
        
        % specific to 2D
        % --------------
        [Yf,Xf] = FreqMesh(fc);
        weights = ( sqrt(2*pi) * sig )^d * exp( -2*pi^2*sig^2 * (Xf.^2 + Yf.^2) );
        % --------------
        
    otherwise
        error('Unknown kernel');
end

end



function [Xg,Yg] = compute_subsampling_grid(fc,model)

switch model.grid
    case 'lattice'
        L = getoptions(model, 'gsize', [64 64]);
        
        q = ceil( (2*fc+1) ./ L );
        t = 1; % TODO!
        
        Lq = L .* q;
        %sgridX  = linspace( 0, (Lq(1)-1)*t, Lq(1) )' / Lq(1);
        sgridX  = linspace( 0, (Lq(1)-1)*t, L(1) )' / Lq(1);
        %sgridY  = linspace( 0, (Lq(2)-1)*t, Lq(2) )' / Lq(2);
        sgridY  = linspace( 0, (Lq(2)-1)*t, L(2) )' / Lq(2);
        [Yg,Xg] = meshgrid(sgridY, sgridX);
        
        %figure, scatter(Xg(:), Yg(:), '.')
        
    case 'boundary'
        L = getoptions(model, 'gsize' , [64, 64]);
        %r = 0.2; %radius TODO: parameter set by user?
        r = getoptions(model, 'gparam', 0.2);
        c = [0.5 0.5]; %center TODO: parameter set by user?
        theta = 2*pi*linspace(0,1,prod(L));
        Xg = c(1) + r*cos(theta);
        Yg = c(2) + r*sin(theta);
        
        
    case 'disc' %TODO
        L = getoptions(model, 'gsize' , [64 64]);
        r = getoptions(model, 'gparam', 0.2    ); % radius
        c = [0.5 0.5]; % center
        
        [Yg,Xg] = ugrid_disc(c,r,L);
        
        %figure, scatter(Xg(:),Yg(:),'.')
        
        
    case 'ring' %TODO
        L = getoptions(model, 'gsize' , [64 64]  );
        r = getoptions(model, 'gparam', [0.2 0.3]);
        c = [0.5 0.5]; % center
        
        [Yg,Xg] = ugrid_ring(c,r(1),r(2),L);
        
end
        
        

end

