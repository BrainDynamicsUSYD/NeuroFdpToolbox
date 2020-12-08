function [R] = GetBurst(R)
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.
% adjust from get_SWR.m function

GammaBurstEvent.spike_sort_range = 6;
GammaBurstEvent.no_std = 2; % 2 std above the mean 1.5 %%% 1.5
GammaBurstEvent.peak_no_std = 2; % 4 std above the mean 3 %%% 1.5
GammaBurstEvent.burst_min_ms = 15;

dt = R.dt;

% burst detection
[no, steps] = size(R.LFP.LFP{1});
spike_hist = R.spike_hist{1};

burst_min_steps = round(GammaBurstEvent.burst_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
GammaBurstEvent.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
GammaBurstEvent.inside_rate = cell(1,no);
GammaBurstEvent.outside_rate = cell(1,no);
GammaBurstEvent.Hz = zeros(1,no);
GammaBurstEvent.is_burst = zeros(size(R.LFP.LFP_broad));

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
            if mean_diff < 0.05 && std_diff < 0.05 %% && mean_diff2 < 0.1 && std_diff2 < 0.1
                GammaBurstEvent.is_burst(i,:) = is_burst_tmp;
                break
            end
        end
        
        % thresholding
        is_burst_tmp = LFP_gamma_hilbert_tmp - hil_mean_baseline_hist(i,r_iter) > GammaBurstEvent.no_std*hil_std_baseline_hist(i,r_iter);
        
        % Discard transient data
        is_burst_tmp(1:transient_steps) = 0;
        
        % persistence_requirement
        is_burst_tmp = persistence_requirement( is_burst_tmp, burst_min_steps );
        
        % peak height requirement
        peak_min = hil_mean_baseline_hist(i,r_iter) + GammaBurstEvent.peak_no_std*hil_std_baseline_hist(i,r_iter);
        [ is_burst_tmp ] = peak_height_requirement( is_burst_tmp, LFP_gamma_hilbert_tmp, peak_min );
        
    end
    
    cutoff = 1;
    [~, burst_du, flat_du, burst_start, ~] = seq_postprocess(GammaBurstEvent.is_burst(i,:), 1, cutoff);
    GammaBurstEvent.Hz(i) = length(burst_du)/(R.dt*R.step_tot*10^-3); % burst frequency in whole simulation
    GammaBurstEvent.burst_du_steps{i} = burst_du;
    GammaBurstEvent.flat_du_steps{i} = flat_du;
    GammaBurstEvent.burst_start_steps{i} = burst_start;
    
    % below four command lines could be deleted
    s_tmp = R.ExplVar.LFP_range_sigma;
    spike_sort_neurons = R.LFP.LFP_neurons{1}(i,:) >= 1/(s_tmp*sqrt(2*pi))*exp(-0.5*(GammaBurstEvent.spike_sort_range/s_tmp)^2);
    
    % firing rate in and out of burst
    GammaBurstEvent.inside_rate{i} = transpose(sum(spike_hist(spike_sort_neurons, is_burst_tmp), 2)/(R.dt*sum(is_burst_tmp)*10^-3));
    GammaBurstEvent.outside_rate{i} = transpose(sum(spike_hist(spike_sort_neurons, ~is_burst_tmp), 2)/(R.dt*sum(~is_burst_tmp)*10^-3));
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
%
% scales = R.LFP.wavelet.scales;
% pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies
% peak.gamma_freq = cell(1,no);
% peak.gamma_freq_step = cell(1,no);
% peak.gamma_wl_amp = cell(1,no);
% peak.gamma_raw_amp = cell(1,no);
% peak.gamma_raw_amp_step = cell(1,no);
% peak.gamma_amp_no_std = cell(1,no);
% for i = 1:no
%     x_tmp = R.LFP.LFP_gamma(i,:);
%     coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
%
%     for j =  1:length(GammaBurstEvent.burst_du_steps{i})
%         a_tmp = GammaBurstEvent.burst_start_steps{i}(j);
%         l_tmp = GammaBurstEvent.burst_du_steps{i}(j);
%
%         [max_coffs_tmp, freq_ind_tmp] = max(coeffs_tmp(a_tmp:a_tmp+l_tmp-1, :), [], 2);
%         [gamma_wl_amp, ind_tmp] = max(max_coffs_tmp);
%         peak.gamma_freq{i} = [peak.gamma_freq{i} pseudoFreq(freq_ind_tmp(ind_tmp))]; % in each burst's peak frequency
%         peak.gamma_freq_step{i} = [peak.gamma_freq_step{i} a_tmp+ind_tmp-1];
%         peak.gamma_wl_amp{i} = [peak.gamma_wl_amp{i} gamma_wl_amp];
%
%         [gamma_raw_amp, gamma_raw_amp_step] = max(R.LFP.LFP_gamma(i, a_tmp:a_tmp+l_tmp-1));
%         peak.gamma_raw_amp{i} = [peak.gamma_raw_amp{i} gamma_raw_amp];
%         peak.gamma_raw_amp_step{i} = [peak.gamma_raw_amp_step{i} a_tmp+gamma_raw_amp_step-1];
%
%     end
%     peak.gamma_amp_no_std{i} =  (peak.gamma_raw_amp{i} - GammaBurstEvent.hil_mean_baseline(i)) / GammaBurstEvent.hil_std_baseline(i);
%
% end
%
% R.LFP.wavelet.peak = peak;
R.LFP.GammaBurstEvent = GammaBurstEvent;
end