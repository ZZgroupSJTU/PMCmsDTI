function [im_mag_spa,im_phase_spa,ParamsUse]=PMC_DTI_spa_recon(ksp,im_init,sens,params)

%% Local low rank constrain to reconstruct the multi-shot diffusion data
% input:  ksp  - [nx,ny,ns,nc,nESP,nshot,ndir,nrep]
%         sens - [nx,ny,ns,nc,nESP]
% Output: im_mag_spa -[nx,ny,ns,1,1,nshot,ndir,nrep]
% Hu Y, Wang X, Tian Q, Yang G, Daniel B, McNab J, et al. Magnetic Resonance in Medicine. 2020;83:1596-1609

%%
pa.lambda1 = 0.8e5; % Regularization parameter on the magnitude images (will be set later based on lambda1).
pa.lambda2 = 0.2; % Regularization parameter on the phase images
pa.winSize = [8,8]; % Block size for the locally low-rank term
pa.phase_update = true; % Whether to update the phase
pa.phase_cycyling = true; % Whether to do phase cycyling (useful when having phase constraint to solve the phase wrapping problem)
pa.iter = 100; % Number of iterations
pa.warmup = true; % Whether to use warm up or not. If used, the first 20% iterations will only have gradient update (no proximal operator update).

%%
if nargin<4
    params=[];
end
[ParamsUse] = UpdateParams(pa, params, true);
%%
[nx,ny,ns,~,~,nshot,ndir,nrep]=size(ksp);
im_mag_spa=zeros(nx,ny,ns,ndir,nrep);
im_phase_spa=zeros(nx,ny,ns,nshot,ndir,nrep);

p = gcp('nocreate');
if ~isempty(p)
    parfor slice_idx=1:ns
        disp(['------------Begin Xinx5 of slice # ',num2str(slice_idx),' ----------------']);tic;
         im_slice=   squeeze(im_init(:,:,slice_idx,:,:,:,:,:,:));
         ksp_slice=  squeeze(ksp(:,:,slice_idx,:,:,:,:,:,:));
         sens_slice= repmat(squeeze(sens(:,:,slice_idx,:)),[1 1 1 ndir]);

         [mag,phase] = Xinx5(im_slice,ksp_slice, sens_slice,...
                                        ParamsUse.iter, ParamsUse.lambda1, ParamsUse.lambda2, ...
                                            ParamsUse.winSize, ParamsUse.phase_update,ParamsUse.phase_cycyling, ParamsUse.warmup);   
                % mag: nx-ny-1-ndir-nex
                % phase: nx-ny-nshot-ndir-nex
                % error: iter-nex     
        im_mag_spa(:,:,slice_idx,:,:,:)=mag;  % [Nx Ny Ns Nb_dir Nrep]
        im_phase_spa(:,:,slice_idx,:,:,:)=phase;
        disp(['------------End Xinx5 of slice # ',num2str(slice_idx),' in ',num2str(toc),' ------------']); 

    end
else
    for slice_idx=1:ns
        disp(['------------Begin Xinx5 of slice # ',num2str(slice_idx),' ----------------']);tic;
         im_slice=   squeeze(im_init(:,:,slice_idx,:,:,:,:,:,:));
         ksp_slice=  squeeze(ksp(:,:,slice_idx,:,:,:,:,:,:));
         sens_slice= repmat(squeeze(sens(:,:,slice_idx,:)),[1 1 1 ndir]);

         [mag,phase] = Xinx5(im_slice,ksp_slice, sens_slice,...
                                        ParamsUse.iter, ParamsUse.lambda1, ParamsUse.lambda2, ...
                                            ParamsUse.winSize, ParamsUse.phase_update,ParamsUse.phase_cycyling, ParamsUse.warmup);   
                % mag: nx-ny-1-ndir-nex
                % phase: nx-ny-nshot-ndir-nex
                % error: iter-nex     
        im_mag_spa(:,:,slice_idx,:,:,:)=mag;  % [Nx Ny Ns Nb_dir Nrep]
        im_phase_spa(:,:,slice_idx,:,:,:)=phase;
        disp(['------------End Xinx5 of slice # ',num2str(slice_idx),' in ',num2str(toc),' ------------']); 

    end
end
