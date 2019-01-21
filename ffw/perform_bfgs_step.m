function [U,exitflag,output] = perform_bfgs_step(U,F,G,options)
%PERFORM_BFGS_STEP bfgs step in ffw algorithm
%   Minimize the non-convex function F(U). G is the gradient wrt U

deal2 = @(varargin) deal(varargin{1:nargout});

nb_cols = size(U,2);

% transform complex variables into real ones
reshc = @(u) reshape( u(1:length(u)/2    ), [length(u)/(2*nb_cols), nb_cols] ) ...
    + 1i *   reshape( u(length(u)/2+1:end), [length(u)/(2*nb_cols), nb_cols] );

flatc = @(Z) [real(Z(:)); imag(Z(:))];

Fb = @(Z) F( reshc(Z) );
Gb = @(Z) flatc( G( reshc(Z) ) );

%figure(10), checkgradient(Fb, Gb, flatc(U)); drawnow;

BFGS_Grad = @(Z) deal2( Fb(Z), Gb(Z) );

[U,fS,exitflag,output] = minFunc(BFGS_Grad, flatc(U), options);
U = reshc(U);


end

