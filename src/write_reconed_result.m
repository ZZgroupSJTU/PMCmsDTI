function write_reconed_result(reconed_result,reconed_dcmF,raw_dcmF,params)
ns = params.dsizebart{6};
ndir = params.dsizebart{7};
nrep = params.dsizebart{8};

% the dicom info extracted from PMCoff msDWI 
[metadata,DicomImg]=Dicom2metadata(raw_dcmF);
dicom_meta=reshape(metadata,[ns,ndir,nrep]);

im_multi_slice_norm = abs(reconed_result(:,:,:))./max(abs(reconed_result(:)));
im_multi_slice_norm = uint16(im_multi_slice_norm*4096);
N_volume = size(im_multi_slice_norm,3);

dicom_meta_multi_slice = permute(dicom_meta(:,:),[2 1]);
im2dicom(im_multi_slice_norm,dicom_meta_multi_slice(:),reconed_dcmF)


%  dcm2niix
dcm_cmd = [' dcm2niix -f %p_%s_%t ', reconed_dcmF];
unix(dcm_cmd);

% save bvals and bvecs 
bvals=params.diffusion.bvals;
bvecs=params.diffusion.bvecs;
valsfile=[reconed_dcmF,filesep,'bvals.txt'];
vecsfile=[reconed_dcmF,filesep,'bvecs.txt'];
save(valsfile, 'bvals', '-ascii'); 
save(vecsfile, 'bvecs', '-ascii'); 

disp('------------End of write_reconed_result------------'); 