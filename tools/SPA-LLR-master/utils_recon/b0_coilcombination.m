function [im_ini] = b0_coilcombination(klb_cc, extraL)
% Coil combination of the fully-sampled b=0 k-space data.
% Input:
%   klb_cc: nx-ny-nc
%   extraL: number of extra ky lines (for partial Fourier)
% Output:
%   im_ini: coil combined image (nx-ny)

acs = extraL*2;
ny = size(klb_cc, 2);
sens = klb_cc.*0;
sens(:,ny/2-acs/2+1:ny/2+acs/2, :) = klb_cc(:,ny/2-acs/2+1:ny/2+acs/2, :); % central k-space
sens = ifft2c(sens); % sensitivity map from the low-resolution image

sens = sens ./ repmat(sos(sens),[1 1 size(sens,3)]); 
im_ini = sum(ifft2c(klb_cc) .* conj(sens),3);
end

