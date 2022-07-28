function [im1,phase] = POCSMUSE(k0,sens,iter,im0,flag)
% This is the implementation of POCSMUSE for multi-shot DWI reconstruction. 
% While there is one thing to be fixed: phase update needs to be included as in POCS-ICE.

% Input:
%   k0: nx-ny-nc-ns
%   sens: nx-ny-1-nc
%   iter: number of iterations
%   im0: initialized image (take the phase as intialization otherwise using SENSE)
%   flag: for intialization whether to use low-resolution phase (0) or not (1).

if nargin < 3
	iter = 200;
end
if nargin < 5
	flag = 0;
end

[p.nx,p.ny,p.nc,p.nshot] = size(k0);

%% Phase estimation
if nargin < 4
    im_step1 = zeros(p.nx,p.ny,p.nshot);
    for ns = 1 : p.nshot
        % phase estimation using SENSE
        im_step1(:,:,ns) = bart('pics -l2 -r 0.01 -w 1',permute(k0(:,:,:,ns),[1 2 4 3]),sens);
    end

    mask = repmat(tri_window(p.nx,p.ny,0),[1 1 p.nshot]); 
    imlow_step1 = ifft2c(fft2c(im_step1).*mask); % low-resolution phase 

else
    mask = repmat(tri_window(p.nx,p.ny,0),[1 1 p.nshot]);
    if flag == 0
	   imlow_step1 = ifft2c(fft2c(im0).*mask);
    else 
	   imlow_step1 = im0;
    end
end

phase = imlow_step1 ./ abs(imlow_step1);

%% Reconstruction. 
% Here we are treating each shot as a coil, and the corresponding sensitivity map is the multiplication of 
% origianl sensitivity map and the phase of that shot. Then we do a POCS-SENSE reconstruction.

sens2 = repmat(squeeze(sens),[1 1 1 p.nshot]).*repmat(permute(phase,[1 2 4 3]),[1 1 p.nc 1]);

im1 = recon_SENSE2D(k0(:,:,:),sens2(:,:,:),iter,1);

end

