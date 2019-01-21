function [UX] = Upsamp1(X,u)
%UPSAMP1 Upsampling operator in 1D
%   UX = UPSAMP1(X,u) upsamples X by a factor u (zeros are inserted)

s = size(X,1);

UX = zeros([u*s, size(X,2)]);
UX(1:u:end,:) = X;

end

