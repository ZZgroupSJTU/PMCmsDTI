function [mag,phase, error,im_llr]=run_spa_llr(ksp,Nshot)
% ksp:  nx-ny-nc-nshot-ndir.
%%
addpath(fullfile('/home/zzgroup/Downloads/bart-0.7.00', 'matlab'));
setenv('TOOLBOX_PATH', '/home/zzgroup/Downloads/bart-0.7.00');
setenv('OMP_NUM_THREADS','1');
%%
senstmp=bart('ecalib -r 32 -c 0.97',permute(ksp(:,:,:,:,1),[1 2 4 3]));
sens=senstmp(:,:,:,:,1);
%%
mask=permute(repmat(eye(4),[size(ksp,2)/Nshot 1]),[3 1 4 2]);
ksp_shot=bsxfun(@times,ksp,mask);
%%
iter = 100; % number of iterations
lambda = 5000; % regularization parameter for LLR term
im_llr=[];
for m=1:size(ksp,5)
%     im_llr(:,:,:,m) = POCSMUSE(ksp_shot(:,:,:,m),sens);
        k_bart=permute(ksp_shot(:,:,:,:,m),[1 2 6 3 5 4]);
        comm = sprintf(['llr = squeeze(bart(',char(39),...
            'pics -R L:7:7:%d -w 1 -i %d',char(39),', k_bart,sens));'], lambda,iter);
%         % Using the pics function in BART to solve the rereconstruction problem. 
%         % Please refer to BART about how to call it. 
        eval(comm);

        im_llr(:,:,:,m)=llr;
end
%%
sens_all=repmat(squeeze(sens),[1 1 1 7]);
%% Setting of reconstruction parameters
% clear all
pa.nslice = 1; % Number of slices
pa.ndir = 6; % Number of diffusion-encoding directions (not including b=0)
pa.nex = 1; 
pa.lambda1 = 0.8e5; % Regularization parameter on the magnitude images (will be set later based on lambda1).
pa.lambda2 = 0.01; % Regularization parameter on the phase images
pa.winSize = [8,8]; % Block size for the locally low-rank term
pa.do_allnex = false; % Whether or not to joint reconstruct all nex
pa.b0 = []; % This is a list storing all the non-diffusion-weighted scans (not including the first b=0 scan). Usually it is empty.
% If using Qiyuan's tensor file for acquisition, set this to "1:15:p.ndir", since the first of every fifteen directions would be the non-diffusion-weighted scan.

pa.b1 = [1:pa.ndir]; % This is a list stroing all diffusion-weighted directions which need to be reconstructed.
pa.b1(pa.b0) = []; 

pa.phase_update = true; % Whether to update the phase
pa.phase_cycyling = true; % Whether to do phase cycyling (useful when having phase constraint to solve the phase wrapping problem)
pa.iter = 100; % Number of iterations
pa.warmup = false; % Whether to use warm up or not. If used, the first 20% iterations will only have gradient update (no proximal operator update).

%pa.do_tv = false; % Apply total variation (1) or L1-wavelet (0) on the phase images. 
%%
[mag,phase, error] = Xinx5(im_llr,...
    ksp_shot, sens_all, pa.iter, pa.lambda1, pa.lambda2, ...
    pa.winSize, pa.phase_update,pa.phase_cycyling, pa.warmup);   
    % mag: nx-ny-1-ndir-nex
    % phase: nx-ny-nshot-ndir-nex
    % error: iter-nex       