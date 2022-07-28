function im_com = recon_SENSE2D(k, smap, iter, flag_plot)
%% reconstruction using SENSE (POCS)
% Input :
%   k(Nx-Ny-Nc): under-sampled kspace 
%   smap : sensitivity map
%   iter : number of iteration
%   flag_plot: whether to plot the data consistency term or not.
% Output : 
%   im_com (Nx Ny): combined image

% By Yuxin Hu, April 22,2015

if (nargin < 3)
    iter = 200;
end
if (nargin < 4)
    flag_plot = 0;
end

mask = k.*0;
mask(abs(k)~=0) = 1;

[Nx Ny Nc] = size(smap);
div = repmat(sum(abs(smap).^2,3),[1 1 Nc]);
div(abs(div) < eps) = 1;
im_iter = ifft2c(k);

for i = 1 : iter
    %disp(i)
    im_iter = im_iter.*conj(smap)./div;
    %im_iter(isnan(im_iter)) = 0; 

    im_com = sum(im_iter,3); % coil combined 
    im_iter = repmat(im_com,[1 1 Nc]);
    k_iter = fft2c(smap.*im_iter);
    k_iter2 = k_iter.*(1-mask)+k; 
    if(flag_plot)
        error2 = abs(k_iter2 - k_iter);
        error(i) = sum(error2(:).^2);
    end
    im_iter = ifft2c(k_iter2);

%     keydown = waitforbuttonpress;
end
if(flag_plot)
    figure,plot(error),title('Reconstruction error','fontsize',14)
    xlabel('number of iteration','fontsize',14)
    ylabel('recon error','fontsize',14)
    im_iter = im_iter.*conj(smap)./div;
    im_iter(isnan(im_iter)) = 0; 
    im_com = sum(im_iter,3);
end

end

