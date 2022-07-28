function newdata = GCC_compress(rawdata,A)
% Compress the k-space data based on the compression matrix.
% Input:
% 	rawdata(kx,ky,kz,nc)

rawdata=fftshift(ifft(fftshift(rawdata,1),[],1),1);

disp('compressing data into virtual coils.');
for slice = 1:size(rawdata,1)
%	disp(slice)
	for j=1:size(A,3) % compressed coil index    
        for k=1:size(rawdata,4) % original coil index
            mat(:,:,k,j)=repmat(A(slice,k,j),[size(rawdata,2) size(rawdata,3) 1]);
        end
        	newdata(slice,:,:,j)=sum(permute(rawdata(slice,:,:,:),[2 3 4 1]).*mat(:,:,:,j),3);
    end
end

% do the circshift back
newdata=fftshift(fft(fftshift(newdata,1),[],1),1);

end

