function [im_allcs,phase] = POCS_ICE(k,sens,iter,lambda,k0,flag)
% POCS-ICE reconstruction of the multi-shot DWI data.
% Input:
%   k : kspace data, nx - ny - nc - ns
%   sens : sensitivity map, nx - ny - nc
%   iter: number of iterations
%   lambda: controls the ratio of orignal and new k-space (between 0 to 1).
%   k0: initialized k-space (the same size as k).
%   flag: whether to use original phase (1) or low-resolution phase (0) in each iteration.
%   
% Output:
%   im_allcs : final output (combined all shots and coils), nx - ny

if nargin < 3 
    iter = 400;
end

if nargin < 4
    lambda = 1;
end
if nargin < 5
    k0 = k.*0;
end
if nargin < 6
    flag = 0;
end

debug = 0;

[nx,ny,nc,ns] = size(k);

div = repmat(sum(abs(sens).^2,3),[1 1 ns]);
div(div == 0) = 1;

mask = k.*0;
mask(abs(k)~=0) = 1; % sampling pattern

mask3 = repmat(tri_window(nx,ny),[1 1 ns]); %  mask for low-resolution phase (triangular window suggested by the POCS-ICE paper)

ktemp = k0;
im_allcs = 0;

for num = 1 : iter
    % Combine all coils for each shot, nx-ny-ns.
    im = ifft2c(ktemp); % image of all coils and all shots, nx-ny-nc-ns.
    im_alls = squeeze(sum(im.*repmat(conj(sens),[1 1 1 ns]),3))./(div+eps); 

    % Calculate phase for each shot 
    if flag == 1 
	    im_alls_low = im_alls;
    else 
        im_alls_low = ifft2c(fft2c(im_alls).*mask3);
    end
    phase = im_alls_low ./ (abs(im_alls_low)+eps);

    % Combine all shots afte removing the phase, nx-ny.
    im_allcs = im_allcs * (1-lambda) + lambda * mean(im_alls.*conj(phase),3); 
    
    % Get the coil-combined image of each shot
    im_alls = repmat(im_allcs,[1 1 ns]) .* phase;

    % Get the kspace of each coil
    im = repmat(permute(im_alls,[1 2 4 3]),[1 1 nc 1]) .* ...
        repmat(sens,[1 1 1 ns]);
    % Update the kspace data
    ktemp = k.*mask + fft2c(im) .* (1-mask);
    
    % Plot the data consistency term.
    if debug
        error2 = abs(fft2c(im).*mask - k);
        error(num) = sum(error2(:).^2);
    end
end

if debug 
    figure,plot(error)
end


