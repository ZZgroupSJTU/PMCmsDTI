function [x_int] = LLe5(llr, do_iniwindow)
%% Intiliaze the input for the SPA-LLR reconstruciton
% Input:
%   llr: initialized complex image, nx-ny-nshot-dir.
% Output:
%   x_int: input for the SPA-LLR reconstruction with size nx-ny-1-(nshot+1)-ndir. The third dimension (leave for coil) is set to 1. The forth dimension is the phase (nshot) and magnitude images (1).

[nx, ny, nshot, ndir] = size(llr);
if nargin < 2
    do_iniwindow = 0;
end

x_int = zeros(nx, ny, (nshot+1), ndir);

if do_iniwindow >= 0
	mask = repmat(tri_window(nx,ny,do_iniwindow),[1 1 nshot ndir]); % Default using hanning window
else
	mask = 1;
end
im_low = ifft2c(mask.*fft2c(llr)); % Low-resolution image
phase = im_low ./ (abs(im_low)+eps); % Low-resolution phase
x_int(:,:,1:nshot,:) = angle(phase);
x_int(:,:,end,:) = abs(mean(llr.*conj(phase),3));




end

