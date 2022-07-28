function [mag, phase,error] = sepMagPhaseRecon3( ...
    x_int,ksp,...
    ESP_op, PFT_op, ...
    MAG_thresh, PHASE_thresh, ...
    MAG_lambda, PHASE_lambda, ...
    MAG_step, PHASE_step, ...
    iter, do_phase_cycling, ...
    do_phase_update, do_warmup)
% Separamete reconstruction of magnitude and phase images for multi-shot multi-direction DWI (SPA-LLR). 
% Input:
%   x_int: initilized input including phase and magnitude, nx-ny-1-(nshot+1)-ndir.
%   ksp: multi-shot multi-coil, multi-direction k-space data, nx-ny-nc-nshot-ndir.
%   ESP_op: sensitivity encoding operator.
%   PFT_op: sampling operator.
%   MAG_thresh: proximal operator of the regularization parameter on the magnitude images.
%   PHASE_thresh: proximal operator of the regularization parameter on the phase images.
%   MAG_lambda: regularization parameter on the magnitude images.
%   PHASE_lambda: regularization parameter on the phase images.
%   MAG_step: step size for magnitude update.
%   PHASE_step: step size for phase update.
%   iter: number of iterations.
%   phase_cycyling: whether to do phase cycyling (useful when having phase constraint to solve the phase wrapping problem)
%   phase_update: whether to update the phase
%   warmup: whether to use warm up or not. If used, the first 20% iterations will only have gradient update (no proximal operator update).

% Output:
%   mag: nx-ny-1-ndir
%   phase: nx-ny-nshot-ndir
%   error: iter

if nargin < 14 
    do_warmup = true;
end


%% Initialize the mag and phase 
mag = x_int(:,:,1,end,:); % nx-ny-1-1-ndir
phase = x_int(:,:,1,1:end-1,:); % nx-ny-1-nshot-ndir
% Again, the 3rd dimension (is 1 now) is for coil dimension.

nshot = size(ksp,4);
ndir = size(ksp,5);
error = zeros(iter,1);

%% Start iteration without regularization (for warm-up)
if do_warmup
for it = 1:iter/5
    % In this for-loop, we are simply doing gradient update for the
    % magnitude image of each direction, phase image of each direction and
    % shot. These steps follow our SPA-LLR paper.
    mag_old = mag;
    phase_old = phase;
    expPp = exp(1j * phase_old);
    alpha_m = MAG_step;
    alpha_p = PHASE_step / (max(abs(mag_old(:)))^2 + eps); % scale the step size to make the update for the magnitude and phase consistent.
    
    % Update residue, mag, and phase
    for d = 1 : ndir
        for n = 1 : nshot
            temp = ksp(:,:,:,n,d) - PFT_op(n) * (ESP_op(d) * (mag_old(:,:,1,1,d).* expPp(:,:,1,n,d)));
            error(it) = error(it) + sum(abs(temp(:)).^2);
            resid(:,:,:,n,d) = ESP_op(d)' * (PFT_op(n)'*temp);
        end
    end
    
    mag = mag_old + alpha_m * real(sum(conj(expPp) .* resid,4));
    
    if do_phase_update
        phase = phase_old + alpha_p * imag(  repmat(mag_old,[1 1 1 nshot]) .* conj(expPp) .* resid  );
    end

end
end

%% Start iteration
for it = (do_warmup*iter/5+1):iter
    % In this for-loop, we are updating the magnitude and phase images
    % based on the gradient and the regularization parameter.
    % These steps follow our SPA-LLR paper.
    
    mag_old = mag;
    phase_old = phase;
    
    expPp = exp(1j * phase_old);
    alpha_m = MAG_step;
    alpha_p = PHASE_step / (max(abs(mag_old(:)))^2 + eps);
    
    % Get random phase wraps
    if (do_phase_cycling)
        r = rand * 2*pi ;
%         phase_wrap = angle( exp(1j * (x_int(:,:,:,1:end-1) + r) ) ) - x_int(:,:,:,1:end-1) - r;
        phase_wrap = angle( exp(1j * (phase_old + r) ) ) - phase_old - r;
    else
        phase_wrap = zeros(size(phase));
    end
    
    % Update residue, mag, and phase based on the gradient.
    for d = 1 : ndir
        for n = 1 : nshot
            temp = ksp(:,:,:,n,d) - PFT_op(n) * (ESP_op(d) * (mag_old(:,:,1,1,d).* expPp(:,:,:,n,d)));
            error(it) = error(it) + sum(abs(temp(:)).^2);
            resid(:,:,:,n,d) = ESP_op(d)' * (PFT_op(n)'*temp);
        end
    end
    
    % Magnitude image update based on the proximal operator
    mag = MAG_thresh(mag_old + alpha_m * real(  sum(conj(expPp) .* resid,4) ), alpha_m * MAG_lambda);
    % All directions are jointly updated.
    
    % Phase image update based on the proximal operator
    if do_phase_update
        for d = 1 : ndir
            for n = 1 : nshot
            phase(:,:,:,n,d) = PHASE_thresh( phase_old(:,:,:,n,d) + alpha_p * imag(mag_old(:,:,1,1,d) .* conj(expPp(:,:,:,n,d)) .* resid(:,:,:,n,d))...
                + phase_wrap(:,:,:,n,d), alpha_p * PHASE_lambda) - phase_wrap(:,:,:,n,d);
            end % loop over shot
        end % loop over direction
    end
   

end


