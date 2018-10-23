function [new_ids] = fft_to_colex(fc)
%FFT_TO_COLEX Switch sorting of moment orders between (*) -fc, ... fc
%(standard colexicographical) and (**) 0,-1, ...,-fc,fc, ... 1 
%(fft-friendly). It can be used in arbitrary dimension.
% 
%   fc:         vector (of size d) of cutoff frequencies in each dimension
%   new_ids:    table (of size (2*fc+1)^d) containing indices such that
%               if moments T are sorted as (**), T(order) is sorted as (*).
%
% CAUTION: as far as I know, it only works if fc(1) = ... = fc(d) -->
% apparently not, from xp, works in any case...

% TODO: remove fc, work only with n

n = 2*fc+1;

new_ids = ifftshift( reshape( (prod(n):-1:1)', [n, 1] ) );
%new_ids = ifftshift( reshape( (1:prod(n))', [n, 1] ) );
% prod? when fc is not the same alogn each dimension???


end

