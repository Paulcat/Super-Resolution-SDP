function [A,As,AsA,Amat,Asmat,Smat] = approximationOperator(fc,model)
%APPROXIMATIONOPERATOR Fourier approximation operator/matrix
%   [A,AS] = APPROXIMATIONOPERATOR(FC,MODEL) returns the approximation
%   operator A and its adjoint AS of the forward operator specified in
%   MODEL, in the Fourier basis [exp(-2i*pi*<k,.>)], -FC <= k <= FC.
%
%   [...,S] = APPROXIMATIONOPERATOR(FC,MODEL) additionally returns the
%   operator S = A'A.
%
%   [...,AMAT,ASMAT] = APPROXIMATIONOPERATOR(FC,MODEL) also returns the 
%   corresponding matrices in the canonical bases (?).
%
%   A takes as input a vector of size PROD(2*FC+1) and returns a vector of
%   size PROD(N), where N depends of the MODEL.
%
%   MODEL has the following fields:
%       FOP     - convolution | subsampled-convolution | foveation
%
%       KERNEL  - type   - Dirichlet | Gaussian
%               - std    - standard deviation (for Gaussian)
%
%       GRID    - shape    - lattice | circle | disc | ring
%               - size     - size of the grid
%               - radius   - radius (for circle)
%               - radii    - radii (for ring)
%
%   For further details, see [Catala, P., Duval, V. and Peyre, G., A
%   low-rank approach to off-the-grid sparse super-resolution]


n = 2*fc+1;
d = length(fc);

% helper
flat = @(x) x(:);
reshn = @(x) reshape(x,n);


[Yf,Xf] = FreqMesh(fc);
%[Xg,Yg] = compute_subsampling_grid(fc,model);
%figure, scatter(Xg(:),Yg(:),'.');

kernel = model.kernel;


switch model.fop
    case 'convolution'
        w = compute_weights(fc,kernel);
        
        A  = @(z) flat( w .* reshn(z) );
        As = @(p) flat( conj(w) .* reshn(p) );
        
        AsA  = @(z) flat( abs(w).^2 .* reshn(z) );
        
        
        if 1
            Amat  = diag(w(:));
            Asmat = diag(w(:))';
            Smat = diag( abs(w(:)).^2 );
        end
        
    case 's-convolution'
        grid   = model.grid;

        w       = compute_weights(fc,kernel);
        [Xg,Yg] = compute_grid(fc,grid);
        L       = getoptions(grid, 'size', [64 64]);
        
        C = 1/ prod(L)^(1/d); % TODO: gamma
        
        if max(2*fc) >= min(L)
            error('grid is insufficiently large wrt fc (TODO)');
        end
        
        if 1
            Amat  = exp( 2i*pi * (Xg(:)*Xf(:)' + Yg(:)*Yf(:)') ) * diag(w(:));
            Asmat = C * Amat';
            Smat  = Asmat * Amat;
        end
        
        if strcmp(grid.shape, 'lattice') % fft implementation
            % helpers
            reshL = @(x) reshape(x, L);
            
            q  = ceil( (2*fc+1) ./ L );
            t  = 1; % TODO!
            Lq = L .* q;
            eL = exp( -2i*pi* ( fc(1)*Xg + fc(2)*Yg ) );
            
            Pad   = @(x) Pad2(fc, x, Lq);
            %Restr = @(x) Restr2(x, [n,1]);
            Restr = @(x) x(1:n(1),1:n(2));
            Ssamp = @(x) Subsamp2(x,t*q);
            Usamp = @(x) Upsamp2(L,x,t*q);
            
            
            %A  = @(z) prod(Lq) * flat( Ssamp( eL .* ifft2( Pad( w.*reshn(z) ) ) ) ); % prod(Lq) counterbalances ifft2
            A = @(z) prod(Lq) * flat( Ssamp( eL .* ifft2( Pad( w(:).*z(:) ) ) ) );
            As = @(p) C * flat( conj(w) .* Restr( fft2( conj(eL) .* Usamp(p(:)) ) ) );
            
            %AsA = @(z) C * flat( abs(w).^2 .* reshn(z) );
            AsA = @(z) 1/C * flat( abs(w).^2 .* reshn(z) ); % TODO: document!!!!!! reason for 1/C is factor N^2 appearing in sum when calculing A'A...
            Smat = C * diag( abs(w(:)).^2 );
            
        else % no fast implementation, store the entire matrix
            %A  = @(z) reshape( Amat  * ( w(:) .* z(:) )    , [L,1] );
            %AS = @(p) reshape( conj(w(:)) .* (ASmat * p(:)), [n,1] );
            
            %A  = @(z) Amat * ( w(:) .* z(:) );
            %As = @(p) conj(w(:)) .* (Asmat * p(:));
            A  = @(z) Amat * z(:);
            As = @(p) Asmat * p(:);
            AsA  = @(z) Smat * z(:);
        end
        
    case 'foveation'
        grid = model.grid;
        
        [Xg,Yg]  = compute_grid(fc,grid);
        L        = getoptions(grid, 'size' , [64 64]);
        cov_type = getoptions(kernel, 'cov', 'default');
        
        if max(2*fc) >= min(L)
            error('grid is insufficiently large wrt fc (TODO)');
        end
        
        Nf = 100; % fft approximation (as Riemann sum)
        [Y,X] = meshgrid( (0:(Nf-1))' / Nf );
        
        switch kernel.mode
            
            case 1 % phig = @(xg,yg) exp( -1/2 * iSig(X,Y,X-sx,Y-sy) )
                iCov = covmat(cov_type);
        
                [J,I] = meshgrid(-1:1);
                I = reshape(I,[1,1,numel(I)]);
                J = reshape(J,[1,1,numel(J)]);
                
                freq1 = (Nf/2+1-fc(1)) : (Nf/2+1+fc(1));
                freq2 = (Nf/2+1-fc(2)) : (Nf/2+1+fc(2));
                
                if max(L) <= 128
                    % full storage: quickly out of memory
                    % P = sum( exp( -1/2 * iSig(X(:)'+I, Y(:)'+J, Xg(:)-X(:)'-I, Yg(:)-Y(:)'-J) ), 3);
                    
                    P = zeros(prod(L),Nf^2);
                    for k=1:numel(I)
                        x = X(:)'+I(k); y = Y(:)'+J(k);
                        E = exp(-1/2 * iCov(x, y, Xg(:)-x, Yg(:)-y));
                        P = P + E;
                    end
                    P = reshape(P,[],Nf,Nf); % matrix [phi(g,x)] + periodization
                    
                    Amat = 1/Nf^2 * fftshift(fftshift(fft(fft(P,[],2),[],3),2),3); % fft2 of each "line" of P
                    Amat = reshape(Amat(:, freq1, freq2 ),[],prod(2*fc+1)); % first (2fc+1) coefficients
                    Amat = conj(Amat);
                    
                else % out of memory issue
                    phij1 = @(xg,yg) sum(exp( -1/2 * iCov(X+I, Y+J, xg-X-I, yg-Y-J) ),3);
                    
                    Amat = zeros(prod(L), prod(2*fc+1));

                    wb = waitbar(0,'Computing matrix A');
                    for j=1:prod(L)
                        uj = phij1( Xg(j), Yg(j) );
                        Fj = 1/Nf^2 * fftshift( fft2(uj) );
                        Fj = Fj(freq1, freq2);

                        Amat(j,:) = conj(Fj(:));
                        
                        %figure(1), surf(X,Y,uj,'LineStyle','none'); view(2); colorbar; drawnow;
                        %figure(1), surf(Xf,Yf,real(Fj),'LineStyle','none'); view(2); drawnow;
                        % normalize?

                        waitbar(j/prod(L),wb);
                    end
                    close(wb);
                end
                
            case 2 % phig = @(xg,yg) exp( -1/2 * iSig(sx,sy,X-sx,Y-sy) )
                [~,Cov,Det] = covmat(cov_type);
                
                c  = 2*pi * sqrt(Det(Xg(:),Yg(:)));
                em = c .* exp(-2i*pi*(Xg(:)*Xf(:)' + Yg(:)*Yf(:)'));
                Amat = em .* exp(-2*pi^2 * Cov(Xg(:),Yg(:),Xf(:)',Yf(:)'));
                Amat = conj(Amat);
                
%                 for j=1:prod(L)
%                     x = Xg(j); y = Yg(j);
%                     
%                     % closed form
%                     em = 2*pi * sqrt(Det(x,y)) * exp(-2i*pi*(x*Xf + y*Yf));
%                     Fj = em .* exp(-2*pi^2 * Sig(x,y,Xf,Yf));
%                     Amat(j,:) = conj(Fj(:));
%                 end
        end
        
%       %APPROXIMATION ERROR????
%         for j=1:prod(L)
%             uj = phij1(Xg(j),Yg(j));
%             Fj = 1/N^2 * fftshift( fft2(uj) );
%             a = sum(sum( abs(Fj).^2 ));
%             Fjc = Fj( (N/2+1-fc) : (N/2+1+fc), (N/2+1-fc) : (N/2+1+fc) );
%             b = sum(sum( abs(Fjc).^2 ));
%             res = res+(a-b);
%         end
        
        Asmat = sqrt(1/prod(L)) * Amat';
        
        %A  = @(z) reshape( Amat  * z(:), [L,1] );
        %AS = @(p) reshape( ASmat * p(:), [n,1] );
        A  = @(z) Amat  * z(:);
        As = @(p) Asmat * p(:);
        
        Smat = Asmat * Amat;
        AsA = @(z) Smat * z(:);

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


function [weights] = compute_weights(fc,kernel)
% COMPUTE_WEIGHTS Fourier coefficients of the kernel

d = numel(fc);

switch kernel.type
    case 'Dirichlet'
        n = 2*fc+1;
        %weights = ones(prod(n), 1) / prod(n); % normalization l1
        weights = ones(n) / prod(n); % normalization l1
        
    case 'Gaussian'
        sig = getoptions(kernel, 'cov', .05);
        
        % specific to 2D
        % --------------
        [Yf,Xf] = FreqMesh(fc);
        weights = ( sqrt(2*pi) * sig )^d * exp( -2*pi^2*sig^2 * (Xf.^2 + Yf.^2) );
        % --------------
        
        %weights = weights / prod(fc);
        
    otherwise
        error('Unknown kernel');
end

end



function [Xg,Yg] = compute_grid(fc,grid)

switch grid.shape
    case 'lattice'
        L = getoptions(grid, 'size', [64 64]);
        
        q = ceil( (2*fc+1) ./ L );
        t = 1; % TODO!
        
        Lq = L .* q;
        %sgridX  = linspace( 0, (Lq(1)-1)*t, Lq(1) )' / Lq(1);
        sgridX  = linspace( 0, (Lq(1)-1)*t, L(1) )' / Lq(1);
        %sgridY  = linspace( 0, (Lq(2)-1)*t, Lq(2) )' / Lq(2);
        sgridY  = linspace( 0, (Lq(2)-1)*t, L(2) )' / Lq(2);
        
        sgridX = linspace(0, L(1)-1, L(1))' / L(1);
        sgridY = linspace(0, L(2)-1, L(2))' / L(2);
        [Yg,Xg] = meshgrid(sgridY, sgridX);
        
        %figure, scatter(Xg(:), Yg(:), '.')
        
    case 'circle'
        L = getoptions(grid, 'size' , [64, 64]);
        %r = 0.2; %radius TODO: parameter set by user?
        r = getoptions(grid, 'radius', 0.2);
        c = [0.5 0.5]; %center TODO: parameter set by user?
        theta = 2*pi*linspace(0,1,prod(L));
        Xg = c(1) + r*cos(theta);
        Yg = c(2) + r*sin(theta);
        
        
    case 'disc' %TODO
        L = getoptions(grid, 'size' , [64 64]);
        r = getoptions(grid, 'radius', 0.2    ); % radius
        c = [0.5 0.5]; % center
        
        [Yg,Xg] = ugrid_disc(c,r,L);
        
        %figure, scatter(Xg(:),Yg(:),'.')
        
        
    case 'ring' %TODO
        L = getoptions(grid, 'size' , [64 64]  );
        r = getoptions(grid, 'radii', [0.2 0.3]);
        c = [0.5 0.5]; % center
        
        [Yg,Xg] = ugrid_ring(c,r(1),r(2),L);
        
end

end

