function [supp,amp,modulus,R] = mvprony(mm,d,options)
%MVPRONY Multivariate Prony extraction
%   [SUPP,AMP] = MVPRONY(MM,D,OPTIONS) computes the support SUPP and
%   amplitudes AMP of a discrete measures matching moments (under moment
%   matrix form) MM, in dimension D.
%
%   SUPP(i,:) gives the location of the i-th source, with amplitude AMP(i)
%
%   Ordering of MM is colexicographical (??) TODO

n = size(mm,1)^(1/d); % n = fc+1 or n = 2fc+1
n = ceil(n); % HACK
n2 = size(mm,1);

options.null = 0;
tol = getoptions(options, 'tol', 1e-3);
shift_mode = getoptions(options, 'shift_mode', 'harmouch');
mode_debug = getoptions(options, 'mode_debug', 0);
smoothing_step = getoptions(options,'smoothing_step',0);

mycellfun = @(fun,cell) cellfun(fun,cell,'UniformOutput',false);

% d-dimensional meshgrid
G = cell(1,d);
[G{:}] = ndgrid(1:n);
G = mycellfun(@(M)M(:), G);

% restriction to avoid out-of-bound
R = mycellfun(@(v)v<n, G);
R = find( prod( cat(2,R{:}), 2 ) );

% re-indexation
G = mycellfun(@(v)v-1, G);
G{1} = G{1}+1; % G contains u, v-1, w-1, etc... meshgrid
nd = n.^(0:d-1);

% compute shifting operators
Shift = cell(1,d);
for i=1:d
    GR = mycellfun(@(v) v(R), G);
    I = sum( nd .* cat(2,GR{:}), 2);
    
    Gi = GR;
    Gi{i} = Gi{i}+1;
    J = sum( nd .* cat(2,Gi{:}), 2);
    
    %Shift{i} = sparse(I,J,ones(length(R),1), n*n, n*n); % shift matrix
    Shift{i} = sparse(I,J,ones(length(R),1), n2, n2); % shift matrix
end

N = cell(1,d); % multiplication matrices
switch shift_mode
    case 'gab'
        % heuristic -- seems to work fine
        [U,S] = myeig(mm,tol);
        for i=1:d
            N{i} = U'*Shift{i}*U;
        end
        
    case 'harmouch'
        % commutes exactly in the sparse case where the hierarchy collapses
        [U,S] = myeig(mm(R,R),tol);
        for i=1:d
            Mi = Shift{i} * mm; Mi = Mi(R,R);
            %N{i} = diag(1./S) * U' * Mi * U;
            N{i} = diag(sqrt(1./S)) * U' * Mi * U * diag(sqrt(1./S));
        end
        
end

% in theory they should commute if hierachy collapses
if mode_debug
    e = norm(N{1}*N{2}-N{2}*N{1}, 'fro') / norm(N{1}*N{2}, 'fro');
    fprintf('Relative commutation error: %.3f\n', e);
end

% joint diagonalization of Nx and Ny
lambda = getoptions(options, 'lambda', [1 -1]);
lambda = lambda(1:d); % HACK!
lambda = reshape(lambda, [1 1 d]);

Nco = sum(lambda .* cat(3, N{:}),3);

% % HACK
% Nco = N{1} - N{2}*N{2};
% %

[H,h] = eig(Nco); h = diag(h);

%norm(Nco*H - H*diag(h),'fro')
%A = (N{1}*H) ./ H;
%B = (N{2}*H) ./ H;


clear options1;
options1 = struct;
switch smoothing_step
    case 'gab'
        options1.H = H;
        options1.tau = 0.09;
        options1.niter = 1000;
        H = joint_diagonalize(N,options1);
        
    case 'cardoso'
        As = cell2mat(N);
        As = reshape(As,[size(Nco),d]);
        H = jeigen_pcg(As);
        H = inv(H);
end

%norm(Nco*H - H*diag(h),'fro')
%A = (N{1}*H) ./ H;
%B = (N{2}*H) ./ H;
%norm(B -diag(diag(B)),'fro')

% eigenvalues: Ni*H = diag(ei)*H
for i=1:d
    %ei = diag(H'*N{i}*H) ./ diag(H'*H);
    ei = diag(inv(H)*N{i}*H);
    %abs(ei)
    
    % map to position on the circe
    %supp(:,i) = mod( -angle(ei)/(2*pi), 1); % TODO: check "minus"
    supp(:,i) = mod( -angle(ei)/(2*pi), 1);
    modulus(:,i) = abs(ei);
end



% compute amplitude: solve Vandermonde system
% TODO: comment
signed = getoptions(options,'signed',1);
if 0
if signed
    ell = [n-1, n-1] / 2; % TODO: remve /2 ??
    %nmom = prod(4*ell+1);
    [Yf,Xf] = meshgrid(-ell(1):ell(1),-ell(2):ell(2));
    Fsupp = exp(-2i*pi* (Xf(:)*supp(:,1)' + Yf(:)*supp(:,2)') );

    [~,~,ids] = marginals(ell,options);
    z = mm(ids);
    z = reshape(z, 4*ell+1);
    z = fftshift(z);

    i0 = 2*ell(1)+1; j0 = 2*ell(2)+1; % TODO: remove?
    z = z(i0-ell(1):i0+ell(1), j0-ell(2):j0+ell(2)); % TODO: remove?
else
    ell = [n-1, n-1]; % TODO handle case ell(1) != ell(2)
    %nmom = prod(ell+1);
    %[Yf,Xf] = meshgrid(0:2*ell(1),0:2*ell(2));
    [Yf,Xf] = meshgrid(-ell(1):ell(1),-ell(2):ell(2));
    Fsupp = exp(-2i*pi * (Xf(:)*supp(:,1)' + Yf(:)*supp(:,2)') );
    
    [~,~,ids] = marginals(ell,options);
    z = mm(ids);
    z = reshape(z,2*ell+1);
    z = fftshift(z);
end


amp = real(Fsupp \ z(:));
end

amp = 0;





clear R;
R.H = H;
R.N = N;
R.h = h;
R.U = U;

end


function [U,S] = myeig(M,tol)

[U,S,V] = eig(M); S = diag(S);
I = find(abs(S)>tol*max(abs(S)));
U = U(:,I); V = V(:,I); S = S(I);

end