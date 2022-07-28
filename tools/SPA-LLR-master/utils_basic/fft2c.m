function res = fft2c(x)

S = size(x);
x = reshape(x,S(1),S(2),prod(S(3:end)));

res = fftc(x,1);
res = fftc(res,2);
res = reshape(res,S);
