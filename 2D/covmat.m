function [Ci,C,D] = covmat(type)
%COVMAT Compute covariance matrix
% SIG = COVMAT(TYPE) returns a covariance matrix SIG of the form
%   SIG = T' * D * T, where
%
%   T = |cos(t) -sin(t)| and D = |d1 0|
%       |sin(t)  cos(t)|         |0 d2|
% t, d1 and d2 may depend on the position (x,y)


switch type
    case 'scalar'
        Ci = @(x,y,v1,v2) (v1.^2 + v2.^2) / .05^2;
        C  = @(x,y,v1,v2) (v1.^2 + v2.^2) * .05^2;
        D  = @(x,y) .05^4;
        
    case 'diagonal'
        r = @(x,y) sqrt( (mod(x,1)-0.5).^2 + (mod(y,1)-0.5).^2 );
        d1 = @(x,y) (.04 + .03 * r(x,y)).^2;
        d2 = @(x,y) (.04)^2;
        %d1 = @(x,y) .05^2;
        %d2 = @(x,y) .05^2;
        
        Ci = @(x,y,v1,v2) v1.^2./d1(x,y) + v2.^2./d2(x,y);
        C  = @(x,y,v1,v2) d1(x,y).*v1.^2 + d2(x,y).*v2.^2;
        D  = @(x,y) 1./d1(x,y)./d2(x,y);
        
    case 'default'
        t = @(x,y) atan2(mod(y,1)-0.5, mod(x,1)-0.5);
        r = @(x,y) sqrt( (mod(x,1)-0.5).^2 + (mod(y,1)-0.5).^2 );
        
        d1 = @(x,y) (.04 + .03 * r(x,y)).^2;
        d2 = @(x,y) (.04)^2;
        
        % inverse covariance matrix
        a = @(x,y) cos(t(x,y)).^2 ./ d1(x,y) + sin(t(x,y)).^2 ./ d2(x,y);
        b = @(x,y) sin(t(x,y)).^2 ./ d1(x,y) + cos(t(x,y)).^2 ./ d2(x,y);
        c = @(x,y) (1./d2(x,y) - 1./d1(x,y)) .* cos(t(x,y)) .* sin(t(x,y));
        
        Ci = @(x,y,v1,v2) a(x,y).*v1.^2 + b(x,y).*v2.^2 + 2*c(x,y).*v1.*v2;
        
        
        % covariance matrix
        d = @(x,y) d1(x,y) .* cos(t(x,y)).^2 + d2(x,y) .* sin(t(x,y)).^2;
        e = @(x,y) d1(x,y) .* sin(t(x,y)).^2 + d2(x,y) .* cos(t(x,y)).^2;
        f = @(x,y) (d2(x,y) - d1(x,y)) .* cos(t(x,y)) .* sin(t(x,y));
        
        C = @(x,y,v1,v2) d(x,y).*v1.^2 + e(x,y).*v2.^2 + 2*f(x,y).*v1.*v2;
        
        
        % determinant
        D = @(x,y) d1(x,y) .*  d2(x,y);
        
    otherwise
        error('TODO')
end


%

% T1 = 100; T2 = 100;
% if 0
%     CC = zeros(2,2,T1,T2);
%     for i=1:T1
%         for j=1:T2
%             CC(:,:,i,j) = [a(i,j), c(i,j); c(i,j), b(i,j)];
%         end
%     end
% end
% 
% if 0
%     CCi = zeros(2,2,T1,T2);
%     for i=1:T1
%         xi = (i-1)/T1;
%         for j=1:T2
%             yj = (j-1)/T2;
%             CCi(:,:,i,j) = [a(xi,yj), c(xi,yj); c(xi,yj), b(xi,yj)];
%         end
%     end
% end
