function A = GCC_A(rawdata, numVrec, calsize,correction)
% Coil compression matrix calculation
% Input:
% 	rawdata (kx,ky,kz,coil)
% 	numVrec: number of virtual coils
% 	calsize: size of the calibration region
%	correction: whether to perform the alignment

if nargin < 3
    calsize = 24;
end
if nargin < 4
    correction = 1;
end

% now go to hybrid space (x,ky,kz,coil)
rawdata=fftshift(ifft(fftshift(rawdata,1),[],1),1);

disp('initialization of coil compression matrices.');
for slice = 1 : size(rawdata,1);
	% find the autocalibration signal
	for i = 1 : size(rawdata,4)	
        tmp = crop(permute(rawdata(slice,:,:,i),[2 3 4 1]),calsize);
		caldata(:, i)=tmp(:);		
	end
	
	[U,S,V]=svd(caldata,'econ');
	% save the SVD results
	svdmatrices(slice,:,:)= V(:,1:numVrec);

end

if(correction)
disp('aligning of coil compression matrices.');
% find out the rotation matrices from the center to both sides
for slice = size(rawdata,1)/2-1:-1:1

	V1 = squeeze(svdmatrices(slice+1,:,:));
	V2 = squeeze(svdmatrices(slice,:,:));

	A = V1'*V2;
	[UA,SA,VA] = svd(A,'econ');
	P = VA*UA';
	svdmatrices(slice,:,:)= (squeeze(svdmatrices(slice,:,:))*P);
end


for slice = size(rawdata,1)/2+1:size(rawdata,1)	
	V1 = squeeze(svdmatrices(slice-1,:,:));
	V2 = squeeze(svdmatrices(slice,:,:));		
	A = V1'*V2;
	[UA,SA,VA] = svd(A,'econ');
	P = VA*UA';
	svdmatrices(slice,:,:)= (squeeze(svdmatrices(slice,:,:))*P);
end
A = svdmatrices;

end

