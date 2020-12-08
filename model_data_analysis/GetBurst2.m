function [R] = GetBurst2(R)
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.
% adjust from GetGamma.m and GetBurst.m functions to get the burst quickly;
% now it is used for generating 2D spatial pattern, i.e. get rid of min
% length limit.
tic;
fprintf('\t Getting Gamma band filtered signal...\n');
dt = 0.1;
fs = 1/(dt*1e-3); % sampling frequency (Hz)

LFP = R.LFP.LFP{1};
[no,steps] = size(R.LFP.LFP{1});
% [no,steps] = size(R.LFP.LFP_gamma_hilbert_abs);

% Butterworth filter
order = 4; % 4th order
lowFreq = 30; % gamma band
hiFreq = 80;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
gaus_width = 12.5; % ms
[ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );
LFP_gamma_hilbert_abs = zeros(no,steps);
for i = 1:no
%     LFP_gamma = LFP(i,:);
    LFP_gamma = filter(b,a,LFP(i,:)); 
    % hilbert transformation & gaussian smoothing
%     LFP_gamma_hilbert_abs(i,:) = abs(hilbert(LFP_gamma));
    LFP_gamma_hilbert_abs(i,:) = conv(abs(hilbert(LFP_gamma)), Kernel,'same'); 
end

% Output results
% R.LFP.LFP_gamma = LFP_gamma;
% R.LFP.gauss_width = gaus_width;
R.LFP.LFP_gamma_hilbert_abs = LFP_gamma_hilbert_abs;
clear LFP_gamma LFP_gamma_hilbert_abs

fprintf('\t Finishing get_Gamma...\n');

GammaBurstEvent.no_std = 2; % 2 std above the mean 1.5
GammaBurstEvent.peak_no_std = 2; % 4 std above the mean 3
% GammaBurstEvent.burst_min_ms = 15;

% burst detection
% burst_min_steps = round(GammaBurstEvent.burst_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
GammaBurstEvent.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
GammaBurstEvent.Hz = zeros(1,no);
GammaBurstEvent.is_burst = zeros(size(R.LFP.LFP_gamma_hilbert_abs));

iter_num = 100;
hil_mean_baseline_hist = zeros(no,iter_num);
hil_std_baseline_hist = zeros(no,iter_num);
for i = 1:no
    LFP_gamma_hilbert_tmp = R.LFP.LFP_gamma_hilbert_abs(i,:);
    is_burst_tmp = zeros(size(LFP_gamma_hilbert_tmp));
    for r_iter = 1:iter_num
        hil_mean_baseline_hist(i,r_iter) =  mean(LFP_gamma_hilbert_tmp(~is_burst_tmp));
        hil_std_baseline_hist(i,r_iter)  =  std(LFP_gamma_hilbert_tmp(~is_burst_tmp));
        
        if r_iter > 1
            mean_diff = abs(hil_mean_baseline_hist(i,r_iter) - hil_mean_baseline_hist(i,r_iter-1))/hil_mean_baseline_hist(i,r_iter-1);
%             mean_diff2 = abs(hil_mean_baseline_hist(i,r_iter) - hil_mean_baseline_hist(i,r_iter-1));
            std_diff = abs(hil_std_baseline_hist(i,r_iter) - hil_std_baseline_hist(i,r_iter-1))/hil_std_baseline_hist(i,r_iter-1);
%             std_diff2 = abs(hil_std_baseline_hist(i,r_iter) - hil_std_baseline_hist(i,r_iter-1));
            if mean_diff < 0.05 && std_diff < 0.05  %% && mean_diff2 < 0.1 && std_diff2 < 0.1
                GammaBurstEvent.is_burst(i,:) = is_burst_tmp;
                break
            end
        end
        
        % thresholding
        is_burst_tmp = LFP_gamma_hilbert_tmp - hil_mean_baseline_hist(i,r_iter) > GammaBurstEvent.no_std*hil_std_baseline_hist(i,r_iter);
        
        % Discard transient data
        is_burst_tmp(1:transient_steps) = 0;
        
        % persistence_requirement
%         is_burst_tmp = persistence_requirement( is_burst_tmp, burst_min_steps );
        
        % peak height requirement
        peak_min = hil_mean_baseline_hist(i,r_iter) + GammaBurstEvent.peak_no_std*hil_std_baseline_hist(i,r_iter);
        [ is_burst_tmp ] = peak_height_requirement( is_burst_tmp, LFP_gamma_hilbert_tmp, peak_min );
        
    end
    
    cutoff = 1;
    [~, burst_du, flat_du, burst_start, ~] = seq_postprocess(GammaBurstEvent.is_burst(i,:), 1, cutoff);
    GammaBurstEvent.Hz(i) = length(burst_du)/(dt*steps*10^-3); % burst frequency in whole simulation
    GammaBurstEvent.burst_du_steps{i} = burst_du;
    GammaBurstEvent.flat_du_steps{i} = flat_du;
    GammaBurstEvent.burst_start_steps{i} = burst_start;
        
end

% GammaBurstEvent.hil_mean_baseline_hist = hil_mean_baseline_hist;
% GammaBurstEvent.hil_std_baseline_hist = hil_std_baseline_hist;
% GammaBurstEvent.hil_mean_baseline = zeros(no,1);
% GammaBurstEvent.hil_std_baseline = zeros(no,1);
% for i = 1:no
%     [~,~,amp_mean_base_tmp] = find(GammaBurstEvent.hil_mean_baseline_hist(i,:), 1, 'last');
%     [~,~,amp_std_base_tmp] = find(GammaBurstEvent.hil_std_baseline_hist(i,:), 1, 'last');
%     GammaBurstEvent.hil_mean_baseline(i) = amp_mean_base_tmp;
%     GammaBurstEvent.hil_std_baseline(i) = amp_std_base_tmp;
% end

R.LFP.GammaBurstEvent = GammaBurstEvent;
toc;
end