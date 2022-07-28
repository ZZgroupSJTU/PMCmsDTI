% This demo loads the b=0 image and the diffusion-weighted images (from
% SPA-LLR), and combine them together into one file for further analysis.
clear all
%% Set acquisition parameters (manually)
ndir = 1 + 6; % Number of total phases (including first b=0, and all diffusion-encoding directions)
b0 = [0]; % Indices for b=0 scans, starting from 0. 
% If using Qiyuan's tensor file for acquisition, set this to "[0, 1:15:ndir-1]", 
% since the first of every fifteen directions would be the non-diffusion-weighted scan.

slice = 1; % Number of slices.

%% Set paths
b1 = [1:ndir];
b1(b0+1) = [];

filepath = '/bmrNAS/people/yuxinh/test_for_github/spa_llr/example_data'; % Path where all data are saved
b0_path = '/prep_sense_gcc'; % Path for b=0 images
md_path = { ...
    '/non_m01p0001i100_sense/', ...
    }; % Path for reconstructed images (from SPA-LLR)


%% Load b=0 images
for d = 1 : length(b0)
    for s = 1 : slice
        if d == 1
            load([filepath,b0_path,'/b0',num2str(s),'.mat']) % from the initial b=0 scan
        else
            load([filepath,b0_path,'/tr',num2str(s + slice * (b0(d)-1)),'.mat']) % from the interleaved b=0 scan
        end            
        im(:,:,s,d) = abs(im_ini);
    end
end
save([filepath,b0_path, '/All_b0.mat'],'im')




%% Load SPA-LLR results
for s = 1 : slice
    load([filepath, md_path{1},'s',num2str(s),'.mat'])
    spa_llr(:,:,s,b1) = abs(mag);
end % loop over all slices
spa_llr(:,:,:,b0+1) = abs(im);

save([filepath,b0_path,'/result'],'spa_llr', 'p', 'pa','md_path','filepath','b0_path')
