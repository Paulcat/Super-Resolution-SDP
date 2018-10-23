function [Xg,Yg] = ugrid_ring(c,r1,r2,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = linspace(0, 1, L(1));
Y = linspace(0, 1, L(2));

R = r1 + r2 * X;
T = 2*pi * Y;

Xg = c(1) + R .* cos(T)';
Yg = c(2) + R .* sin(T)';


%figure, scatter(Xg(:), Yg(:), '.');

end

