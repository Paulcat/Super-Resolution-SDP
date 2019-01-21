function [F,FS] = FourierOperator(fc)
%FOURIEROPERATOR Returns the Fourier operator F and its adjoint FS (Fourier
%synthesis)
%   fc is the cutoff frequency

Xf = FreqMesh(fc);

F  = @(x)    exp( -2i*pi * Xf(:)*x(:)' );
FS = @(DX,p) real( F(DX)' * p(:) );




end

