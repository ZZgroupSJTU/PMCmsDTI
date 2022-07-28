function im2dicom(im,metadata,dicomname)
im=uint16(im(:,:,:));
metadata=metadata(:);
if (length(metadata)<size(im,3))
    error('The size of metadata is smaller than the image');
end

if nargin<3
    dicomname='test';
end

folderName=dicomname;
if (isfolder(folderName))
   rmdir(folderName,'s');
end
mkdir(dicomname);

%%
dicomdict('factory') 
dictfile='/zzNAS/UMR790SJTU/PMCproject/reconstruction/utils/UIH_dicom_dict_backup.txt';
dicomdict('set', dictfile);

for m=1:size(im,3)
    dicomwrite(im(:,:,m), [dicomname,filesep,num2str(m,'%08d'),'.dcm'], metadata{m},'WritePrivate',true, 'CreateMode', 'copy');
end