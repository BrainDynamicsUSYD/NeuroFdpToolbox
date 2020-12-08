function GammaBurstEvent = GetBurstFR(FR)
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.
% adjust from GetBurst.m function
% applying baseline + 2s.d. method on FR

GammaBurstEvent.no_std = 2; % 2 std above the mean 1.5 %%% 1.5
GammaBurstEvent.peak_no_std = 2; % 4 std above the mean 3 %%% 1.5
GammaBurstEvent.burst_min_ms = 0;

dt = 1;

% burst detection
[no, steps] = size(FR);

burst_min_steps = round(GammaBurstEvent.burst_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
GammaBurstEvent.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
GammaBurstEvent.Hz = zeros(1,no);
GammaBurstEvent.is_burst = zeros(no,steps);

iter_num = 100;
hil_mean_baseline_hist = zeros(no,iter_num);
hil_std_baseline_hist = zeros(no,iter_num);
for i = 1:no
    is_burst_tmp = zeros(size(FR));
    for r_iter = 1:iter_num
        hil_mean_baseline_hist(i,r_iter) =  mean(FR(~is_burst_tmp));
        hil_std_baseline_hist(i,r_iter)  =  std(FR(~is_burst_tmp));
        
        if r_iter > 1
            mean_diff = abs(hil_mean_baseline_hist(i,r_iter) - hil_mean_baseline_hist(i,r_iter-1))/hil_mean_baseline_hist(i,r_iter-1);
            std_diff = abs(hil_std_baseline_hist(i,r_iter) - hil_std_baseline_hist(i,r_iter-1))/hil_std_baseline_hist(i,r_iter-1);
            if mean_diff < 0.05 && std_diff < 0.05 
                GammaBurstEvent.is_burst(i,:) = is_burst_tmp;
                break
            end
        end
        
        % thresholding
        is_burst_tmp = FR - hil_mean_baseline_hist(i,r_iter) > GammaBurstEvent.no_std*hil_std_baseline_hist(i,r_iter);
        
        % Discard transient data
        is_burst_tmp(1:transient_steps) = 0;
        
        % persistence_requirement
        is_burst_tmp = persistence_requirement( is_burst_tmp, burst_min_steps );
        
        % peak height requirement
        peak_min = hil_mean_baseline_hist(i,r_iter) + GammaBurstEvent.peak_no_std*hil_std_baseline_hist(i,r_iter);
        [ is_burst_tmp ] = peak_height_requirement( is_burst_tmp, FR, peak_min );
        
    end
    
    cutoff = 1;
    [~, burst_du, flat_du, burst_start, ~] = seq_postprocess(GammaBurstEvent.is_burst(i,:), 1, cutoff);
    GammaBurstEvent.Hz(i) = length(burst_du)/(dt*steps*10^-3); % burst frequency in whole simulation
    GammaBurstEvent.burst_du_steps{i} = burst_du;
    GammaBurstEvent.flat_du_steps{i} = flat_du;
    GammaBurstEvent.burst_start_steps{i} = burst_start;
end

GammaBurstEvent.hil_mean_baseline_hist = hil_mean_baseline_hist;
GammaBurstEvent.hil_std_baseline_hist = hil_std_baseline_hist;
GammaBurstEvent.hil_mean_baseline = zeros(no,1);
GammaBurstEvent.hil_std_baseline = zeros(no,1);
for i = 1:no
    [~,~,amp_mean_base_tmp] = find(GammaBurstEvent.hil_mean_baseline_hist(i,:), 1, 'last');
    [~,~,amp_std_base_tmp] = find(GammaBurstEvent.hil_std_baseline_hist(i,:), 1, 'last');
    GammaBurstEvent.hil_mean_baseline(i) = amp_mean_base_tmp;
    GammaBurstEvent.hil_std_baseline(i) = amp_std_base_tmp;
end
end