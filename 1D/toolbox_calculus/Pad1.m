function [PX] = Pad1(X,PadSize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

PX = zeros(PadSize);
PX(1:size(X,1), :) = X;


end

