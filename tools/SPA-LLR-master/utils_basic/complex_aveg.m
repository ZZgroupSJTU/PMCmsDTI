function im_c = complex_aveg( im )
% Complex average of multiple images 
% im: nx-ny-nex
im = squeeze(im);
[nx,ny,nex] = size(im);
if nex >= 1 
    mask = repmat(tri_window(nx,ny,0),[1 1 nex]);
    im_low = ifft2c(fft2c(im).*mask);
    phase = im_low ./ abs(im_low);
    % removing low-resolution phase then average different nex
    im_c = mean(im.*conj(phase),3);
else
    im_c = im;
end
    
end

