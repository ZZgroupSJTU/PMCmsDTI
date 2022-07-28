function [ktemp] = preProcess(k_hb,slice,r,nshot,scale)
% Given the k-space data of all slices loaded from Orchestra, this function preprocess the data for the following reconstruction.
% The preprocessing includes: (1) split the data into multiple shots with zero-filling (since all shots are combined into one matrix);
% (2) extract the required slice; (3) under-samples the matrix if an acceleration factor is given; (4) scales the data if a scaling factor is given.

% Input: 
% 	k_hb:nx-ny-nc-nex-slice
% 	slice: slice index
% 	r:reduction factor
% 	nshot: number of shots
% 	scale: scaling factor, if none, then use 1
% Output:
%   ktemp: nx-ny-nc-nshot-nex

if nargin < 5
    scale = 1;
end

[nx,ny,nc,nex,~] = size(k_hb);
ktemp = zeros(nx,ny,nc,nshot,nex);

for n = 1 : nex
    for ns = 1 : nshot
        ktemp(:,ns:nshot:ny,:,ns,n) = k_hb(:,ns:nshot:ny,:,n,slice)/scale; 
    end
end


ktemp = ktemp(:,:,:,1:r:end,:);

end

