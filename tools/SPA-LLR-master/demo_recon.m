addpath(genpath(pwd))
%% Setting of reconstruction parameters
% clear all
pa.nslice = 1; % Number of slices
pa.ndir = 6; % Number of diffusion-encoding directions (not including b=0)
pa.nex = 1; 
pa.lambda1 = 0.1; % Regularization parameter on the magnitude images (will be set later based on lambda1).
pa.lambda2 = 0.001; % Regularization parameter on the phase images
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
pa.filepath = 'example_data'; % The path where the initialized SENSE data are saved.

md_path = { ...
    '/non_m01p0001i100_sense', ...
    };
lambda1 = [0.1]; % A list of regularization parameters (on magnitude images) which are going to be used. 
% This makes the grid search of the regularization parameter easier.

for d = 1 : length(md_path)
    if(~exist([pa.filepath, md_path{d}]))
        mkdir([pa.filepath, md_path{d}]);
    end
end
%%
for s = 1 : pa.nslice
    pa.filepath
    fprintf([' reconstructing slice ' num2str(s) ' \n'])
    clear ksp im_llr im_llr2
    for d = 1 : length(pa.b1)
        load([pa.filepath,'/prep_sense_gcc/tr',num2str(s+(pa.b1(d)-1)*pa.nslice),'.mat']);
        ksp(:,:,:,:,d,:) = khb_cc; % nx-ny-nc-nshot-ndir-nex
        im_llr(:,:,:,d,:) = im_ini; % nx-ny-nshot-ndir-nex
        sens_all(:,:,:,d) = sens; % nx-ny-nc-ndir. Different directions may have different sensitivity maps.
    end
    
    if pa.do_allnex 
        % If reconstructing all nex together, then we can simply treat
        % different nex as different shots in one nex. This is achieved by
        % permutation and dimension combination.
        im_llr = permute(im_llr,[1 2 4 3 5]);
        im_llr = permute(im_llr(:,:,:,:),[1 2 4 3]);
        ksp = permute(ksp,[1 2 3 5 4 6]);
        ksp = permute(ksp(:,:,:,:,:),[1 2 3 5 4]);
    end
    
    for d = 1 : length(md_path)
        for n = 1 : pa.nex 

            pa.lambda1 = lambda1(d);
            [mag(:,:,:,:,n),phase(:,:,:,:,n), error(:,n)] = Xinx5(im_llr(:,:,:,:,n),...
                ksp(:,:,:,:,:,n), sens_all, pa.iter, pa.lambda1, pa.lambda2, ...
                pa.winSize, pa.phase_update,pa.phase_cycyling, pa.warmup);   
                % mag: nx-ny-1-ndir-nex
                % phase: nx-ny-nshot-ndir-nex
                % error: iter-nex       
        save([pa.filepath,md_path{d}, '/s',num2str(s),'.mat'],'mag','pa','error','phase');    
        end 
    end % For different regularization parameters
end