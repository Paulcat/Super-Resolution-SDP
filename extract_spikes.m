function [x,z] = extract_spikes(U,fc)
%EXTRACT_SPIKES returns the discrete measure whose moment matrix is UU* 
%   Support extraction is implemented following Algorithm 4.2 in "Moments,
%   positive polynomials and their applications", J.B. Lasserre, 2010.

d = length(fc); % dimension


% * Prune columns of U *
% * ****************** *
%U(:,1)  = []; % 1st column is only zeros
[A,S,B] = svd(U,'econ');
S2      = diag(S).^2;
% r       = sum( S2 / max(S2) > 1e-3); --> too demanding
r       = sum( S2 / max(S2) > 5e-2 );
U = A(:,1:r) * S(1:r,1:r);
%U       = A(:,1:r) * S(1:r,1:r) * B(1:r,1:r); % TODO: not sure truncating
% columns of B is the best way: why the last ones? seems some other choices
% work too for Lasserre extraction, 

%U       = U(:,1:r);

W = rref(U.').'; % column echelon form


% find pivot line-indices
% there are no columns with only zeros
[~,pivot_ids] = max(W~=0,[],1);
r = length(pivot_ids); % TODO: r = nombre de colonnes dans U quoi


pows      = gencolex(fc);
gen_basis = pows(pivot_ids,:);


% Compute multiplication matrices
N = zeros(r,r,d);
for i=1:d
    ei = zeros(size(gen_basis)); ei(:,i) = 1;
    xi_gen_basis = gen_basis + ei;
    [~,pos] = ismember(xi_gen_basis, pows, 'rows');
    %Ni = W(pos,:);
    
    N(:,:,i) = W(pos,:); % ith multiplication matrix (multiplication by xi)
end

% Form random linear combination of multiplication matrices
rweights = rand(d,1); % random weights
rweights = reshape(rweights / norm(rweights,1), [1 1 d]); % summing up to 1
NN = sum( N .* rweights, 3);


% Atoms are located in common eigenspaces
[Q,~] = schur(NN); % Q: is this optimal?
z = zeros(r,d); % --> TODO: plot pour remplacer la figure sur les racines du polynome dual
for i=1:d
    z(:,i) = diag(Q' * N(:,:,i) * Q);
end


% % display z
% t = linspace(0,1,1024);
% for i=1:d
%     figure; hold on;
%     plot(exp(2i*pi*t),'k','LineWidth',1.5);
%     plot(real(z(:,i)), imag(z(:,i)), 'b.', 'MarkerSize', 20);
%     axis equal
%     drawnow;
%     hold off;
% end

% 
theta = 1/2/pi * angle(z);
x = mod( 1 - theta, 1); %




end

