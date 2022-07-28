%% In order to use the bart toolbox
addpath(fullfile('~/bart-0.7.00', 'matlab'));
setenv('TOOLBOX_PATH', '~/bart-0.7.00');
setenv('OMP_NUM_THREADS','1');

%%
addpath(genpath('../../PMC_msDWI'));

%% 
% Data Reconstruction from k-space to image space
% 
%%%%%%%%%%%%%%%
%% Load multi-shot data
load('../Data/ms_dwi/fig5_ms_dwi.mat')

% Dims - [Nro Npe Nslice Ncoils ~ Nshot Nbdirection Nrepetition]
[nx,ny,ns,nc,~,nshot,ndir,nrep]=size(PMCoff.ksp);

%% PMC OFF / LLR OFF
% Combine multi-shot k-space to the full k-space
ksp_shot_combine_pmcoff = sum(PMCoff.ksp,6);

% Transform k-space to image space and coil combination
im_PMCoff_LLRoff = sos(cifftn(ksp_shot_combine_pmcoff,[1 2]),4); 
figure(); imshowMRI(abs(squeeze(im_PMCoff_LLRoff(:,:,16,:,:,:,2:33))),0,[4 8])

%% PMC OFF / LLR ON
% SPA-LLR reconstruction
% Hu Y, Wang X, Tian Q, Yang G, Daniel B, McNab J, et al. Magnetic Resonance in Medicine. 2020;83:1596-1609

% Generate the initial im_llr
im_llr=PMC_DTI_LLR_recon(PMCoff.ksp,PMCoff.sens,PMCoff.params.recon,0);

% SPA-LLR recon
[im_mag_spa,im_phase_spa]=PMC_DTI_spa_recon(PMCoff.ksp,im_llr,PMCoff.sens);
im_PMCoff_LLRon = im_mag_spa;
disp('-------------End of PMC_DTI_spa_LLR_recon-------------------');

figure(); imshowMRI(abs(squeeze(im_PMCoff_LLRon(:,:,16,2:33))),0,[4 8])

%% PMC ON / LLR OFF
% Combine multi-shot k-space to the full k-space
ksp_shot_combine_pmcon = sum(PMCon.ksp,6);

% Transform k-space to image space and coil combination
im_PMCon_LLRoff = sos(cifftn(ksp_shot_combine_pmcon,[1 2]),4); 
figure(); imshowMRI(abs(squeeze(im_PMCon_LLRoff(:,:,16,:,:,:,2:33))),0,[4 8])

%% PMC ON / LLR ON
% SPA-LLR reconstruction
% Hu Y, Wang X, Tian Q, Yang G, Daniel B, McNab J, et al. Magnetic Resonance in Medicine. 2020;83:1596-1609

% Generate the initial im_llr
im_llr=PMC_DTI_LLR_recon(PMCon.ksp,PMCon.sens,PMCon.params.recon,0);

% SPA-LLR recon
[im_mag_spa,im_phase_spa]=PMC_DTI_spa_recon(PMCon.ksp,im_llr,PMCon.sens);
im_PMCon_LLRon = im_mag_spa;
disp('-------------End of PMC_DTI_spa_LLR_recon-------------------');

figure(); imshowMRI(abs(squeeze(im_PMCon_LLRon(:,:,16,2:33))),0,[4 8])

%% 
% Write Reconstructed result to dicom format and convert them to nii files
% 
%%%%%%%%%%%%%%%
%% Dicom -  PMC OFF / LLR ON
% the LLR recon result of PMCoff msDW
im_pre2dcm_PMCoff = permute(im_PMCoff_LLRon,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder of raw PMCoff msDWI from scanner
dcmF_PMCoff_LLRoff='../Data/ms_dwi/epi_dti_pmc_4s_step_off_dk_7401';

% the dicom folder for saving the recon result of PMCoff msDWI 
dcmF_PMCoff_LLRon = '../Data/ms_dwi/PMCoff_msDWI_ReconedDcm';

% writing dicom
write_reconed_result(im_pre2dcm_PMCoff,dcmF_PMCoff_LLRon,dcmF_PMCoff_LLRoff,PMCoff.params);

%% Dicom - PMC ON / LLR ON
% the LLR recon result of PMCon msDW
im_pre2dcm_PMCon = permute(im_PMCon_LLRon,[1 2 4 3 5]);

% the dicom folder of raw PMCon msDWI from scanner
dcmF_PMCon_LLRoff='../Data/ms_dwi/epi_dti_pmc_4s_step_on_5001';

% the dicom folder for saving the recon result of PMCon msDWI 
dcmF_PMCon_LLRon = '../Data/ms_dwi/PMCon_msDWI_ReconedDcm';

% writing dicom
write_reconed_result(im_pre2dcm_PMCon,dcmF_PMCon_LLRon,dcmF_PMCon_LLRoff,PMCon.params);

%% 
% DTI analysis by FSL Toolbox
% Require FSL installed and GPU
%%%%%%%%%%%%%%%
%% FSL path setting 
setenv( 'FSLDIR', '/usr/local/fsl');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); 
clear fsldir fsldirmpath;
%% FSL DTI -  PMC OFF / LLR ON
mk_dicom_folder_off=dcmF_PMCoff_LLRon;
volume_total_num = PMCoff.params.dsizebart{7}*PMCoff.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_off,volume_total_num);

%% FSL DTI -  PMC ON / LLR ON
mk_dicom_folder_on=dcmF_PMCon_LLRon;
volume_total_num = PMCon.params.dsizebart{7}*PMCon.params.dsizebart{8};% ndir*nrep

fsl_dtifit(mk_dicom_folder_on,volume_total_num);




