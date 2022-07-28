function [ mag,phase,error ] = Xinx5(im_llr,ktemp_cc,sens,iter,lambda1,lambda2,winSize,do_phase_update,do_phase_cycling, do_warmup,do_iniwindow)
%% Multi-direction non-convex reconstruction

% Input:
%   im_llr: initilized complex images, nx-ny-nshot-ndir.
%   ktemp_cc: nx-ny-nc-nshot-ndir.
%   sens: nx-ny-nc-ndir.
%   iter: number of iterations.
%   lambda1: regularization parameter on the magnitude images.
%   lambda2: regularization parameter on the phase images.
%   winSize: block size for the locally low-rank term
%   phase_update: whether to update the phase
%   phase_cycyling: whether to do phase cycyling (useful when having phase constraint to solve the phase wrapping problem)
%   warmup: whether to use warm up or not. If used, the first 20% iterations will only have gradient update (no proximal operator update).
%   do_iniwindow: Window for the initialization images (negative: all-pass, 0: hanning window, 1: triangular window).

% Output:
%   mag: nx-ny-1-ndir
%   phase: nx-ny-nshot-ndir
%   error: iter
%% Default parameters
if nargin < 8
	do_phase_update = true;
end
if nargin < 9
	do_phase_cycling = true;
end
if nargin < 10
	do_warmup = true;
end
if nargin < 11
	do_iniwindow = 0;
end

%% Extract parameters
nshot = size(im_llr,3);
nx = size(im_llr,1);
ny = size(im_llr,2);
nc = size(ktemp_cc,3);
ndir = size(ktemp_cc,5);

%% Intilization
x_int1 = permute(LLe5(im_llr, do_iniwindow),[1 2 5 3 4]); % nx-ny-1-(nshot+1)-ndir 
% The third dimension (leave for coil) is set to 1.
% The forth dimension is the phase (nshot) and magnitude images (1).

%% Set reconstruction parameters

MAGstep = 1*0.9/nshot; % Step size for magnitude update (notice that it is scaled by nshot). 
% Since there is only one magnitude image and nshots phase images.

PHASEstep = 1*0.9; % Step size for phase update.

mask = ktemp_cc ~= 0; % Sampling mask.

%% Create linear operators for the data consistency term

% Sensitivity encoding operator
ESP_op = ESPIRiT(sens);
if size(sens,4) < ndir
	sens = repmat(sens,[1 1 1 ndir]);
end

for n = 1 : ndir
	ESP_op(n) = ESPIRiT(sens(:,:,:,n));
end

% Sampling operator
PFT_op = p2DFT(mask(:,:,:,1,1),[nx, ny, nc]);

% Currently assuming the sampling pattern is the same for all different directions
for n = 1 : nshot
    PFT_op(n) = p2DFT(mask(:,:,:,n,1),[nx, ny, nc]);
end

%% Construct the proximal operator

% MAG_thresh = @(x, lambda) wave_thresh(x, lambda, 'db6');
MAG_thresh = @(x, lambda) LLR_thresh(x, lambda, winSize); % Locally low-rank for magnitude images
% MAG_thresh = @(x, lambda) wave_thresh(x, lambda, 'db8');

PHASE_thresh = @(x, lambda) wave_thresh(x, lambda, 'db8');

%% SPA-LLR reconstruction
[mag, phase, error] = sepMagPhaseRecon3(x_int1, ktemp_cc, ESP_op, PFT_op, ...
                                MAG_thresh, PHASE_thresh,...
                                lambda1, lambda2/nshot, ... % Notice that the regularization parameter for phase images is scaled by a factor of nshot.
                                MAGstep, PHASEstep, ...
                                iter, do_phase_cycling,...
			                    do_phase_update, do_warmup);

end

