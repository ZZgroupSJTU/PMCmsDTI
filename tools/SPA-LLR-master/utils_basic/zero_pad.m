function y = zero_pad(x, dim)
% Zero-pad the matrix. x is the matrix to be padded. And dim is a list which describes the size after zero-padding.
% This function is based on Prof. Miki Lustig's ESPIRiT util code.
% Also notice that we only pad the first two dimensions, unless the matrix is a 5D matrix.

if length(dim) == 4
    [kx ky nc nc] = size(x);
    kx_1 = floor(dim(1)/2-kx/2) + 1;
    kx_2 = kx_1 + kx - 1;
    ky_1 = floor(dim(2)/2-ky/2) + 1;
    ky_2 = ky_1 + ky - 1;
    y = zeros(dim);
    y(kx_1:kx_2,ky_1:ky_2,:,:) = x;
end

if length(dim) == 3
    [kx ky nc] = size(x);
    kx_1 = floor(dim(1)/2-kx/2) + 1;
    kx_2 = kx_1 + kx - 1;
    ky_1 = floor(dim(2)/2-ky/2) + 1;
    ky_2 = ky_1 + ky - 1;
    y = zeros(dim);
    y(kx_1:kx_2,ky_1:ky_2,:) = x;
end

if length(dim) == 5
    [kx ky kz nc nc] = size(x);
    kx_1 = floor(dim(1)/2-kx/2) + 1;
    kx_2 = kx_1 + kx - 1;
    ky_1 = floor(dim(2)/2-ky/2) + 1;
    ky_2 = ky_1 + ky - 1;
    kz_1 = floor(dim(3)/2-kz/2) + 1;
    kz_2 = kz_1 + kz - 1;
    y = zeros(dim);
    y(kx_1:kx_2,ky_1:ky_2,kz_1:kz_2,:,:) = x;
end

end


