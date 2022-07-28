%%
addpath(genpath('../../PMC_msDWI'));

%% Load data
load('../Data/ss_dwi/no_motion/Ref.mat')

% Dims - [Nro Npe Nslice Ncoils ~ Nshot Nbdirection Nrepetition]
[nx,ny,~,nc,~,ndir,nrep]=size(Ref.ksp);

%% Transform k-space to image space
im_raw=flip(imrotate(squeeze(sos(cifftn(Ref.ksp,[1 2]),4)),90),2); % RO PE Ns  NbRep Nrep
figure(); imshowMRI(abs(squeeze(im_raw(:,:,16,2:33))),0,[4 8])

%% 
% Write Reconstructed result to dicom format and convert them to nii files
% 
%%%%%%%%%%%%%%%
%% 1. Dicom -  No motion - Ref
% the recon result of no motion
im_pre2dcm_ref = permute(im_raw,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder from scanner
dcmF_raw='../Data/ss_dwi/no_motion/epi_dti_pmc_1s_nm_on_9101';

% the dicom folder for saving the recon result 
dcmF_recon = '../Data/ss_dwi/no_motion/reconed_dcm';

% writing dicom
write_reconed_result(im_pre2dcm_ref,dcmF_recon,dcmF_raw,Ref.params);

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
%% FSL DTI -  No motion - Ref
mk_dicom_folder_ref=dcmF_recon;
volume_total_num = Ref.params.dsizebart{7}*Ref.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_ref,volume_total_num);

%% 2. For continous motion and PMC ON data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('../Data/ss_dwi/cm_PMCon/cm_PMCon.mat')
[nx,ny,~,nc,~,ndir,nrep]=size(cm_PMCon.ksp);

% Transform k-space to image space
im_raw=flip(imrotate(squeeze(sos(cifftn(cm_PMCon.ksp,[1 2]),4)),90),2); % RO PE Ns  NbRep Nrep
figure(); imshowMRI(abs(squeeze(im_raw(:,:,16,2:33))),0,[4 8])

% the recon result 
im_pre2dcm_cm_PMCon = permute(im_raw,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder from scanner
dcmF_raw='../Data/ss_dwi/cm_PMCon/epi_dti_pmc_1s_cm_on_8701';

% the dicom folder for saving the recon result 
dcmF_recon = '../Data/ss_dwi/cm_PMCon/reconed_dcm';

% writing dicom
write_reconed_result(im_pre2dcm_cm_PMCon,dcmF_recon,dcmF_raw,cm_PMCon.params);

% FSL DTI -  cm_PMCon
mk_dicom_folder_ref=dcmF_recon;
volume_total_num = cm_PMCon.params.dsizebart{7}*cm_PMCon.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_ref,volume_total_num);
%% 3. For continous motion and PMC OFF data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('../Data/ss_dwi/cm_PMCoff/cm_PMCoff.mat')
[nx,ny,~,nc,~,ndir,nrep]=size(cm_PMCoff.ksp);

% Transform k-space to image space
im_raw=flip(imrotate(squeeze(sos(cifftn(cm_PMCoff.ksp,[1 2]),4)),90),2); % RO PE Ns  NbRep Nrep
figure(); imshowMRI(abs(squeeze(im_raw(:,:,16,2:33))),0,[4 8])

% the recon result 
im_pre2dcm_cm_PMCoff = permute(im_raw,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder from scanner
dcmF_raw='../Data/ss_dwi/cm_PMCoff/epi_dti_pmc_1s_cm_off_8901';

% the dicom folder for saving the recon result 
dcmF_recon = '../Data/ss_dwi/cm_PMCoff/reconed_dcm';

% writing dicom
write_reconed_result(im_pre2dcm_cm_PMCoff,dcmF_recon,dcmF_raw,cm_PMCoff.params);

% FSL DTI -  cm_PMCon
mk_dicom_folder_cm_PMCoff=dcmF_recon;
volume_total_num = cm_PMCoff.params.dsizebart{7}*cm_PMCoff.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_cm_PMCoff,volume_total_num);

%% 4. For step motion and PMC ON data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('../Data/ss_dwi/step_PMCon/step_PMCon.mat')
[nx,ny,~,nc,~,ndir,nrep]=size(step_PMCon.ksp);

% Transform k-space to image space
im_raw=flip(imrotate(squeeze(sos(cifftn(step_PMCon.ksp,[1 2]),4)),90),2); % RO PE Ns  NbRep Nrep
figure(); imshowMRI(abs(squeeze(im_raw(:,:,16,2:33))),0,[4 8])

% the recon result 
im_pre2dcm_step_PMCon = permute(im_raw,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder from scanner
dcmF_raw='../Data/ss_dwi/step_PMCon/epi_dti_pmc_1s_step_on_8301';

% the dicom folder for saving the recon result 
dcmF_recon = '../Data/ss_dwi/step_PMCon/reconed_dcm';

% writing dicom
write_reconed_result(im_pre2dcm_step_PMCon,dcmF_recon,dcmF_raw,step_PMCon.params);

% FSL DTI -  cm_PMCon
mk_dicom_folder_step_PMCon=dcmF_recon;
volume_total_num = step_PMCon.params.dsizebart{7}*step_PMCon.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_step_PMCon,volume_total_num);
%% 5. For step motion and PMC OFF data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('../Data/ss_dwi/step_PMCoff/step_PMCoff.mat')
[nx,ny,~,nc,~,ndir,nrep]=size(step_PMCoff.ksp);

% Transform k-space to image space
im_raw=flip(imrotate(squeeze(sos(cifftn(step_PMCoff.ksp,[1 2]),4)),90),2); % RO PE Ns  NbRep Nrep
figure(); imshowMRI(abs(squeeze(im_raw(:,:,16,2:33))),0,[4 8])

% the recon result 
im_pre2dcm_step_PMCoff = permute(im_raw,[1 2 4 3 5]); % [Nx Ny Nb_dir Ns*Nrep]


% the dicom folder from scanner
dcmF_raw='../Data/ss_dwi/step_PMCoff/epi_dti_pmc_1s_step_off_8501';

% the dicom folder for saving the recon result 
dcmF_recon = '../Data/ss_dwi/step_PMCoff/reconed_dcm';

% writing dicom
write_reconed_result(im_pre2dcm_step_PMCoff,dcmF_recon,dcmF_raw,step_PMCoff.params);

% FSL DTI -  cm_PMCon
mk_dicom_folder_step_PMCoff=dcmF_recon;
volume_total_num = step_PMCoff.params.dsizebart{7}*step_PMCoff.params.dsizebart{8}; % ndir*nrep

fsl_dtifit(mk_dicom_folder_step_PMCoff,volume_total_num);











