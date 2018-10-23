function [Xg,Yg] = ugrid_disc(c,r,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = linspace(0, 1, L(1));
Y = linspace(0, 1, L(2));
[Yg,Xg] = meshgrid(Y,X);

points     = [Xg(:) Yg(:)];
new_points = zeros(size(points));

%figure;

for i=1:size(points,1)
    x = points(i,1);
    y = points(i,2);
    
    radius = r;
    
    a = 2*x - 1;
    b = 2*y - 1;
    
    if a^2 > b^2
        radius = radius * a;
        phi    = pi/4 * b/a;
    else
        radius = radius * b;
        phi    = pi/2 - pi/4 * a/b;
    end
    
    if (a==0 && b==0)
        new_points(i,1) = c(1);
        new_points(i,2) = c(2);
    else
        new_points(i,1) = c(1) + radius * cos(phi);
        new_points(i,2) = c(2) + radius * sin(phi);
    end

    %scatter(new_points(:,1),new_points(:,2),'.');
    %drawnow;
end

%figure, scatter(new_points(:,1), new_points(:,2), '.');
Xg = reshape( new_points(:,1), L);
Yg = reshape( new_points(:,2), L);

% R = r    * X;
% T = 2*pi * Y;
% 
% Xg = c(1) + R .* cos(T)';
% Yg = c(2) + R .* sin(T)';

end

