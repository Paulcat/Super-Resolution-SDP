function [F,FS] = FourierOperator(fc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Yf,Xf] = FreqMesh(fc);

F = @(x,y) exp( -2i*pi*( Xf(:)*x(:)' + Yf(:)*y(:)' ) );

Y = @(DX,DY) meshgrid(DY,DX);
X = @(DX,DY) meshgrid(DX,DY)';
FS = @(DX,DY,p) numel(DX)*numel(DY) * ...
    real( exp( -2i*pi*( fc(1)*X(DX,DY) + fc(2)*Y(DX,DY) ) ) .* ...
        ifft2( Pad2(p, [numel(DX), numel(DY), 1]) ) );


end

