% Mutli-direction DWI reconstruction: for (1) initialization of SPA-LLR and (2) prep of deep learning reconstruction
% Using GCC for coil compression


%% Load b=0 k-space data (for future sensitivity map calculation).
if (~isempty(p.b0) && p.dir > 15) 
    % Using 15 is because Qiyuan's tensor file interleaves one b0 in every 15 directions.
    file_to_read = p.dir - mod(p.dir-1,15);
    load([p.filename,'/k' int2str(file_to_read+1) '.mat'])
    % If this is used, we will use the neighboor b=0 data for sensitivity map calculation instead of the first one at the very beginning.
else 
    load([p.filename,'/k1.mat'])
end

% k0: nx-ny-nc-nex-slice
if size(k0,4) == 2
    k0(:,:,:,2:2:end,:) = -k0(:,:,:,2:2:end,:); % Interesting here. It seems the second nex is negative of the first nex. Pay attention to this when you want to average of multiple nex.
end
% k_lb = mean(k0,4);

%k_lb = mean(k0(:,:,:,1,:),4); % We simply use the first average here.
k_lb = k0;

%% Extract information from the acqusition header
p.extraL = p0.hnover; % Number of extral lines
p.nshot = p0.ileaves; % Number of shots
p.tSlice = p0.nslices/p0.nphases; % Number of slices
p.partialf = 1; % Partial Fourier is always used
p.ACS = 2*p.extraL; % Number of ACS lines which is double of the extra acquired ky lines.

%% Set reconstruction parameters
p.slice = [1:p.tSlice]; % The slice index to be reconstructed (usually all slices).
p.r = [1]; % Reduction factors (this is for retrospective undersampling, not used).

p.iter = 200; % Number of iterations for shot-LLR reconstruction.

p.allnex = 0; % Whether to reconstruct all nex of each slice and direction together (0 or 1).

p.caverage = false; % Complex average (removing the low-resolution phase) or not after the shot-LLR reconstruction.
% This option will not change much here since we are doing homodyne first, which make the image real-valued.


%% Load diffusion-weighted data
load([p.filename,'/k' int2str(p.dir+1) '.mat']) 
if size(k0,4) == 2
    k0(:,:,:,2:2:end,:) = -k0(:,:,:,2:2:end,:); % Same here. It seems the second nex is negative of the first nex. Pay attention to this when you want to average of multiple nex.
end

k_hb = k0;
[p.nx2,p.ny2,p.nc,p.nex,~] = size(k0); % p.nx2 and p.ny2 save the original matrix size if zero-filled is used.
if exist('NX')
    % Zero-fill the k-space if the size is given.
    k_lb = zero_pad(k_lb,[NX,NY,p.nc,size(k_lb,4),p.tSlice]);
    k_hb = zero_pad(k_hb,[NX,NY,p.nc,p.nex,p.tSlice]); 
end

if p.dir == 1
    % Calculate the scaling factor based on the central slice if this is for the first direction.
    p.scale = max(max(sos(ifft2c(squeeze(k_lb(:,:,:,1,p.slice(round(length(p.slice)/2))))))));
end

k_lb = k_lb / p.scale;
k_hb = k_hb / p.scale;

[p.nx,p.ny,p.nc,p.nex,~] = size(k_hb); % This is the actual matrix size we are going to deal with.

figure('vis','off')

%% Reconstruction of each slice and direction

for s = 1 : length(p.slice)
for r = 1 : length(p.r)
    %% Process the data to be reconstructed
    ktemp = preProcess(k_hb,p.slice(s),p.r(r),p.nshot,1); % nx-ny-nc-nshot-nex
    
    %% Coil compression
    if p.gcc
        [klb_cc,khb_cc] = GCC(squeeze(k_lb(:,:,:,:,p.slice(s))),ktemp, p.v,[40,1]);
        %   klb_cc: compressed b=0 kspace (nx-ny-vc-nex)
        %   khb_cc: compressed diffusion-weighted kspace (nx-ny-vc-nshot-nex)
    else
        klb_cc = squeeze(k_lb(:,:,:,:,p.slice(s)));
        khb_cc = ktemp;
        p.v = p.nc;
    end

    %% Reconstruction
    image_index = (p.dir - 1) * length(p.slice) + s;
    % This is the index to count how many images have been reconstructed. Since we are saving one file for each slice and direction. This number is going to be in the file name.

    % Sensitivity map calculation
    [calib ~] = bart('ecalib -r 20', permute(klb_cc(:,:,:,1),[1 2 4 3])); % Use the first nex for sensitivity map calculation
    sens = squeeze(bart('slice 4 0', calib));
    
    % Deal with the data from b=0 scan
    if p.dir == 1  % only process for the first time
        % k2 = sum(klb_cc, 4);
        k2 = klb_cc;
        for n = 1 : size(k2, 5)
            im_ini(:,:,n) = squeeze(b0_coilcombination(k2(:,:,:,n), p.extraL)); % Coil compression of the b=0 data
            eval(['im_ini(:,:,n) = bart(' char(39) 'homodyne -C -I 1 ' num2str(0.5+p.extraL/p.ny) char(39) ',im_ini(:,:,n));']); % Homodyne reconstruction
        end
        im_ini = complex_aveg(im_ini);
        save([p.savepath, '/b0',num2str(image_index),'.mat'],'k2','sens','im_ini','p')        
    end
    
    % Deal with all the following data
    if ismember(p.dir, p.b0)        
        k2 = sum(khb_cc,4);
        for n = 1 : size(k2, 5)
            im_ini(:,:,n) = squeeze(b0_coilcombination(k2(:,:,:,:,n), p.extraL)); % Coil compression of the b=0 data
            eval(['im_ini(:,:,n) = bart(' char(39) 'homodyne -C -I 1 ' num2str(0.5+p.extraL/p.ny) char(39) ',im_ini(:,:,n));']); % Homodyne reconstruction
        end
        im_ini = complex_aveg(im_ini);
        save([p.savepath, '/tr',num2str(image_index),'.mat'],'k2','sens','im_ini','p')
    else
        if strcmp(p.reconmethod,'LLR')
            % Shot-LLR reconstruction
            [im_ini,LLR(:,:,s,r)] = MSrecon3(khb_cc,sens,p.iter,'CLEAR',0.5+p.extraL/p.ny,p.caverage,p.allnex,p.lambda);
            save([p.savepath, '/tr',num2str(image_index),'.mat'],'khb_cc','sens','im_ini','p')

        elseif strcmp(p.reconmethod,'SENSE')
            for nex = 1 : p.nex
                for ns = 1 : p.nshot
                    % SENSE reconstruction for each shot
                    im_ini(:,:,ns,nex) = bart('pics -l2 -r 0.001 -w 1',permute(khb_cc(:,:,:,ns,nex),[1 2 4 3]),permute(sens,[1 2 4 3]));
                end
            end
            save([p.savepath, '/tr',num2str(image_index),'.mat'],'khb_cc','sens','im_ini','p')
% 
%     else
%         [im_ini,MUSE(:,:,s,r)] = MSrecon3(khb_cc,sens,p.iter,'MUSE',0.5+p.extraL/p.ny,p.caverage,p.allnex,p.lambda);
%         save([p.savepath, '/tr',num2str(image_index),'.mat'],'khb_cc','sens','im_ini','p')

        end
    end
    
end % end of loop over all R
end % end of loop over all slices

% save([p.savepath, '/LLR',num2str(file_index),'_',num2str(p.tSlice),'.mat'],'LLR')
