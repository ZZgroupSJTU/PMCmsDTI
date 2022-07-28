function res = ifft2c(x)

S = size(x);

x = reshape(x,S(1),S(2),prod(S(3:end)));


res = ifftc(x,1);
res = ifftc(res,2);

res = reshape(res,S);

