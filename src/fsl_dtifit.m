function fsl_dtifit(mk_dicom_folder,volume_total_num)

    tmp1 = dir([mk_dicom_folder,filesep,'*.nii']).name;
%     tmp1 = dir([mk_dicom_folder,filesep,'*.nii.gz']).name;
    raw_nii = [mk_dicom_folder,filesep,tmp1];
    %%
    bet_cmd = ['bet ' raw_nii ' ' mk_dicom_folder '/raw_brain -F -f 0.3 -g 0 -m '];
    unix(bet_cmd);
    %%
     tmp11 = dir([mk_dicom_folder,filesep,'*brain.nii.gz']).name;
     raw_brain_nii = [mk_dicom_folder,filesep,tmp11];
%% Eddy
    pre_raw=['cd ' mk_dicom_folder ';mkdir eddy;cd eddy;touch index.txt; touch acqparams.txt; mkdir eddy_result; cd ..'];
    unix(pre_raw);
    addpath(genpath(pwd));    

    % generate mask, since eddy will work in a non-distorted space we will base the mask on my_data.nii.gz
    A1=['fslmaths  ' raw_brain_nii ' -Tmean ' mk_dicom_folder '/eddy/mean_data'];
    unix(A1);
    A3 =['bet ' mk_dicom_folder '/eddy/mean_data ' mk_dicom_folder '/eddy/mean_data_brain -A -f 0.2 -g 0 -m'];
    unix(A3);

    % acqparams.txt
    B=['printf "0 1 0 0.099" > ' mk_dicom_folder '/eddy/acqparams.txt'];
    unix(B);

    % index.txt
    C=['indx=" "; for ((i=1; i<=' num2str(volume_total_num) '; i+=1)); do indx="$indx 1"; done;  echo $indx > ' mk_dicom_folder '/eddy/index.txt']; % volume number = 22
    unix(C);
    
    % required in eddy
    tmp2 = dir([mk_dicom_folder,filesep,'*.bvec']).name;
    bvec = [mk_dicom_folder,filesep,tmp2];
    tmp3 = dir([mk_dicom_folder,filesep,'*.bval']).name;
    bval = [mk_dicom_folder,filesep,tmp3];
    
    D1 =['eddy_cuda9.1 --imain=' raw_brain_nii ' --mask=' mk_dicom_folder '/eddy/mean_data_brain_mask --acqp=' mk_dicom_folder '/eddy/acqparams.txt --index=' mk_dicom_folder '/eddy/index.txt --bvecs=' bvec ' --bvals=' bval ' --repol --cnr_maps --residuals --out=' mk_dicom_folder '/eddy/eddy_result/eddy_corrected_data'];
    unix(D1);
    
    % edge, -out    data_processed/eddy/eddyresult/eddy_brain_inskull_mask.nii.gz  
    E1=['fslmaths  ' mk_dicom_folder '/eddy/eddy_result/eddy_corrected_data -Tmean ' mk_dicom_folder '/eddy/eddy_result/mean_eddy'];
    unix(E1);
    E = ['bet ' mk_dicom_folder '/eddy/eddy_result/mean_eddy ' mk_dicom_folder '/eddy/eddy_result/mean_eddy_brain -A -f 0.2 -g 0'];
    unix(E);
    
%% DTIFIT
    % before eddy
    pre_raw=['cd ' mk_dicom_folder ';mkdir dtifit_before_eddy']; unix(pre_raw);addpath(genpath(pwd));
    F = ['dtifit --data=' raw_brain_nii ' --mask=' mk_dicom_folder '/eddy/mean_data_brain_mask --bvecs=' bvec ' --bvals=' bval ' --out=' mk_dicom_folder '/dtifit_before_eddy/dti_before_eddy'];
    unix(F);
    
    % after eddy
    pre_eddy=['cd ' mk_dicom_folder ';mkdir dtifit_after_eddy']; unix(pre_eddy);addpath(genpath(pwd));
    eddy_nii = [mk_dicom_folder '/eddy/eddy_result/eddy_corrected_data'];
    % generate mask based on eddy_corrected_data.nii.gz
    G1=['fslmaths  ' eddy_nii ' -Tmean ' mk_dicom_folder '/dtifit_after_eddy/mean_eddy_data'];
    unix(G1);

    G3 =['bet ' mk_dicom_folder '/dtifit_after_eddy/mean_eddy_data ' mk_dicom_folder '/dtifit_after_eddy/mean_eddy_data_brain -A -f 0.2 -g 0 -m'];
    unix(G3);

    H = ['dtifit --data=' eddy_nii ' --mask=' mk_dicom_folder '/dtifit_after_eddy/mean_eddy_data_brain_mask --bvecs=' mk_dicom_folder '/eddy/eddy_result/eddy_corrected_data.eddy_rotated_bvecs --bvals=' bval ' --out=' mk_dicom_folder '/dtifit_after_eddy/dti_after_eddy'];
    unix(H);
    disp('------------End of fsl_dtifit.m------------'); 
end