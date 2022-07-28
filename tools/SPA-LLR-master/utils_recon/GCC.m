function [klb_cc,khb_cc] = GCC(k_lb,k_hb,vc,calsize)
% Performs coil compression to one slice of the non-diffusion-weighted and diffusion-weighted data. The compression matrix is calculated based on the 
% fully sampled non-diffusion-weighted data. Geometric coil compression is used. And this code is mainly based on the author's (Dr. Tao Zhang) code.
%
% Input:
%   k_lb: b=0 k-space data (nx-ny-nc-nex).
%   k_hb: diffuion-weighted k-space data (nx-ny-nc-nshot-nex).
%   vc: number of virtual coils to be compressed.
%   calsize: calculation size for GCC, defauly size [40,1].

% output:
%   klb_cc: compressed b=0 kspace (nx-ny-vc-nex).
%   khb_cc: compressed diffusion-weighted kspace (nx-ny-vc-nshot-nex).

if nargin < 4
    calsize = [40,1];
end

A = GCC_A(permute(k_lb(:,:,:,1),[1 2 4 3]), vc, calsize,1); % calculated the compression matrix based on the non-diffusion-weighted data.
for nex = 1 : size(k_lb,4)
    % Compressed the diffusion-weighted data for each slice and nex.
    klb_cc(:,:,:,nex) = squeeze(GCC_compress(permute(k_lb(:,:,:,nex),[1 2 4 3]),A)); % Compressed the b=0 data.
end


khb_cc = k_hb .*0;
khb_cc(:,:,vc+1:end,:,:) = [];
for ns = 1 : size(k_hb,4)
    for nex = 1 : size(k_hb,5)
    	% Compressed the diffusion-weighted data for each slice and nex.
        khb_cc(:,:,:,ns,nex) = squeeze(GCC_compress(permute(k_hb(:,:,:,ns,nex),[1 2 4 3]),A));
    end
end

end

