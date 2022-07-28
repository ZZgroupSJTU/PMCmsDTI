function im_llr=PMC_DTI_LLR_recon(ksp,sens,params,debug)
%% Local low rank constrain to reconstruct the multi-shot diffusion data
% input:  ksp  - [nx,ny,ns,nc,nESP,nshot,ndir,nrep]
%         sens - [nx,ny,ns,nc,nESP]
% Output: im_llr -[nx,ny,ns,1,1,nshot,ndir,nrep]
% Hu Y, Wang X, Tian Q, Yang G, Daniel B, McNab J, et al. Magnetic Resonance in Medicine. 2020;83:1596-1609

disp('-----------------Start LLR reconstruction ----------------------');
tic;
if isfield(params,'lambda')
    lambda=params.lambda;
else
    lambda=1000;
end

if isfield(params,'iter')
    iter=params.iter;
else
    iter=200;
end

if isfield(params,'blksize')
    blksize=params.blksize;
else
    blksize=[7 7];
end

[nx,ny,ns,~,~,nshot,ndir,nrep]=size(ksp);
im_llr=zeros(nx,ny,ns,1,1,nshot,ndir,nrep);
p = gcp('nocreate');
if ~isempty(p)
    parfor slice_idx=1:ns
        for bdir_idx=1:ndir
            for rep_idx=1:nrep
                bartcmd=sprintf('pics -R L:%d:%d:%d -w 1 -i %d', blksize(1),blksize(2),lambda,iter);
                llr=bart(bartcmd,ksp(:,:,slice_idx,:,:,:,bdir_idx,rep_idx),sens(:,:,slice_idx,:));
                im_llr(:,:,slice_idx,1,1,:,bdir_idx,rep_idx)=llr;
            end
        end
    end
else
    for slice_idx=1:ns
        for bdir_idx=1:ndir
            for rep_idx=1:nrep
                bartcmd=sprintf('pics -R L:%d:%d:%d -w 1 -i %d', blksize(1),blksize(2),lambda,iter);
                llr=bart(bartcmd,ksp(:,:,slice_idx,:,:,:,bdir_idx,rep_idx),sens(:,:,slice_idx,:));
                im_llr(:,:,slice_idx,1,1,:,bdir_idx,rep_idx)=llr;
            end
        end
    end
end

if (debug)
    figure();imshow3(abs(im_llr(:,:,round(end/2),1,1,:,round(end/2),1)));
    figure();imshow3(angle(im_llr(:,:,round(end/2),1,1,:,round(end/2),1)));  
end
disp(['LLR reconstruction done in: ', num2str(toc)]);