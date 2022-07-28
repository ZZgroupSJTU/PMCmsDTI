function [imc,im] = MSrecon3(k,sens,iter,method,partialF,aveg,allnex,lambda,k0,l1_flag)
% Reconstruction of multishot DWI with different options. 
% Input: 
%   k(nx-ny-nc-nshot-nex): undersampled kspace
%   sens(nx-ny-nc): sensitivity map. for CLEAR-VCS, its size should be nx-ny-nc-2.
%   iter: number of iterations
%   method: MUSE, ICE, LLR, LLR-VCS, LLR_TWOC. Notice that LLR_TWOC (low rank + sparse) requires our modified BART.
%   partialF: This is the partial Fourier factor. If partial Fourier is not used, then this should be 0.
%   aveg: average methods (0 for real-valued or 1 for complex average).
%   allnex: whether (1) or not (0) recon all nex together.
%   lambda: regularization parameter for LLR.
%   k0: initialization for POCS-MUSE or POCS-ICE. Notice POCS-MUSE asks for image initialization, while POCS-ICE ask for k-space.
%   l1_flag: regularization parameter for the additional L1 contraint in LLR reconstruction (0 for none).

% Output:
%   imc: before homodyne (nx-ny)
%   im: after homodyne (nx-ny)

if nargin < 9
    k0 = k; % the input itself as initialization of POCS-ICE.
end 

if nargin < 10
    l1_flag = 0; % No additional L1-constraint in default.
end

if allnex
    % If reconstructing all nex together, then treat different nex just as differnt shots.
    k = k(:,:,:,:);
    k0 = k0(:,:,:,:);
end

[nx,ny,nc,nshot,nex] = size(k);

switch method
    case 'MUSE'
        im_muse = zeros(nx,ny,nex);
        for n = 1 : nex
            if nargin < 9
                % whether or not to use initlization
                [im_muse(:,:,n),~] = POCSMUSE(k(:,:,:,:,n),permute(sens,[1 2 4 3]),iter);
            else
                [im_muse(:,:,n),~] = POCSMUSE(k(:,:,:,:,n),permute(sens,[1 2 4 3]),iter,k0(:,:,:,n));
            end
        end
        imc = im_muse;
    case 'ICE'
        im_ice = zeros(nx,ny,nex);
        for n = 1 : nex
            [im_ice(:,:,n),~] = POCS_ICE(k(:,:,:,:,n),sens,iter,1, k0(:,:,:,:,n));
        end
        imc = im_ice;

    case 'LLR'
        im_llr = zeros(nx,ny,nshot,nex);
        sens = permute(sens,[1 2 4 3]);
        for n = 1 : nex
            ktemp = reshape(k(:,:,:,:,n),[nx ny 1 nc 1 1 1 1 1 1 nshot]);
	        if l1_flag == 0
                comm = sprintf(['im_llr(:,:,:,n) = squeeze(bart(',char(39),...
                    'pics -R L:7:7:%d -w 1 -i %d',char(39),', ktemp,sens));'], lambda,iter);
	        else 
                comm = sprintf(['im_llr(:,:,:,n) = squeeze(bart(',char(39),...
                    'pics -l1 -r %d -R L:7:7:%d -w 1 -i %d',char(39),', ktemp,sens));'], l1_flag, lambda, iter);
	        end
            eval(comm);
        end
        imc = im_llr(:,:,:,:);
    case 'LLR-VCS'
        im_llr = zeros(nx,ny,2*nshot,nex);

        sens2 = repmat(sens(:,:,:,1),[1 1 1 nshot]);
        sens2(:,:,:,nshot+1:2*nshot) = repmat(sens(:,:,:,2),[1 1 1 nshot]);
        sens2 = reshape(sens2,[nx ny 1 nc 1 1 1 1 1 1 2*nshot]); % different shots/frames have different maps
        
        for n = 1 : nex
            ktemp = k(:,:,:,:,n);
            ktemp(:,:,:,end+1:2*end) = conj(flip(flip(ktemp,1),2));
            ktemp = reshape(ktemp,size(sens2));
            if l1_flag == 0
	            comm = sprintf(['im_llr(:,:,:,n) = squeeze(bart(',char(39),...
                    'pics -R L:7:7:%d -w 1 -i %d',char(39),', ktemp,sens2));'], lambda,iter);
            else 
    
       	        comm = sprintf(['im_llr(:,:,:,n) = squeeze(bart(',char(39),...
                    'pics -l1 -r 0.001 -R L:7:7:%d -w 1 -i %d',char(39),', ktemp,sens2));'], lambda,iter);eval(comm);
	        end
            eval(comm);
        end
        imc = im_llr(:,:,:,:);
    otherwise 
        warning('unexpected method')
end

im = imc(:,:,:);

if (partialF ~= 0) && ~strcmp(method,'CLEAR-VCS')
    for n = 1 : size(im,3)
        eval(['im(:,:,n) = bart(' char(39) 'homodyne -C -I 1 ' num2str(partialF) char(39) ',im(:,:,n));']);
    end
end

if (aveg)
    im = complex_aveg(im);
else
    im = mean(abs(im), 3);
end


end

