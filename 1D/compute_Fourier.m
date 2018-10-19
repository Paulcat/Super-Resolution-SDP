function [F,FS] = compute_Fourier(fc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Xf = FreqMesh(fc);

F  = @(x)    exp( -2i*pi * Xf(:)*x(:)' );
FS = @(DX,p) real( F(DX)' * p(:) );




end

