function [metadata,DicomImg]=Dicom2metadata(DicomPath)

dirOutput=dir(fullfile(DicomPath,'*.dcm'));
DicomFile={dirOutput.name}';
index=1;
DicomImg=[];
for k=1:length(DicomFile)
    DicomInfoOut = dicominfo(fullfile(DicomPath, DicomFile{k})) ;
    if DicomInfoOut.FileSize>100
        metadata{index}=DicomInfoOut;
        DicomImg(:,:,index) = dicomread(fullfile(DicomPath, DicomFile{index})) ;
        index=index+1;
    end
end