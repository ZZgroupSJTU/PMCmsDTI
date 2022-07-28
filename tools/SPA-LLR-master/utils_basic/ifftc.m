function res = ifftc(x,dim)

res = ifftshift(x,dim);
res = ifft(res,[],dim);
res = sqrt(size(x,dim))*fftshift(res,dim);

