function res = fftc(x,dim)

res = ifftshift(x,dim);
res = fft(res,[],dim);
res = 1/sqrt(size(x,dim))*fftshift(res,dim);

