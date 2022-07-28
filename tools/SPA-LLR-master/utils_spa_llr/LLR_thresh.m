function res = LLR_thresh(im, lambda, winSize)
% Proximal operator of the locally low-rank regularization.
% Input:
%   im: nx-ny-nex/frame
%   lambda: regularization parameter
%   winSize: block size for LLR



% Zero-pad the image to a multiple of winSize
im = squeeze(im);
[nx ny nc] = size(im);
bx = ceil(nx/winSize(1));
by = ceil(ny/winSize(2));
zpadx = bx*winSize(1);
zpady = by*winSize(2);

imzpad = zeros(zpadx, zpady, nc);
imzpad(1:nx, 1:ny, :) = im;

% Random shift to avoid the blocky artifacts
nxrs = randperm(winSize(1));
nyrs = randperm(winSize(2));
imzpad = circshift(imzpad,[nxrs(1)-1, nyrs(1)-1, 0]);

% Soft-thresholding of the singular values of each local matrix
for xx = 1 : bx
    for yy = 1 : by
        % Construct the matrix
        imb = imzpad((xx-1)*winSize(1)+1:xx*winSize(1),(yy-1)*winSize(2)+1:yy*winSize(2),:);
        imb = reshape(imb,[winSize(1)*winSize(2),nc]);

        % SVD then soft-thresholding
        [Ub,Sb,Vb] = svd(imb,'econ');
        sdiagb = diag(Sb);
%         sdiagb(1)
        %mu = lambda*sdiagb(1);
        mu = lambda;
        sdiagb = sdiagb - mu;
        sdiagb(sdiagb < 0) = 0;
  
        % Recover the image
        imb = Ub*(diag(sdiagb))*Vb';
        imb = reshape(imb,[winSize(1),winSize(2),nc]);
        imzpad((xx-1)*winSize(1)+1:xx*winSize(1),(yy-1)*winSize(2)+1:yy*winSize(2),:,:) = imb;
    end
end

% Undo the random shift
imzpad = circshift(imzpad,[-nxrs(1)+1, -nyrs(1)+1, 0]);
res = imzpad(1:nx, 1:ny, :);
res = permute(res,[1 2 4 5 3]);
end
