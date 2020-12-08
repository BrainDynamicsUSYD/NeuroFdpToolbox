function GammaBurstEvent = find_Burst_1D(sigIn,fsTemporal,min_time,badChannels,...
    numChannels,stdVal)
% GammaBurstEvent = find_Burst_1D(sigIn,fsTemporal,badChannels, ...
% numberChannels) ;
%
% This function implements the procdure to find burst in 1D data
% 
% Input:
% sigIn          (Channel X Time) signals
% fsTemporal     sample frequency
% numChannels    number of channels, default = 1
%
% Output:
% GammaBurstEvent        structure including:
%    .is_burst           boolean with same size of sigIn for bursts 
%    .burst_du_steps     time duration for bursts
%    .burst_start_steps  start time for bursts
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings
if ~exist('badChannels', 'var') 
    nanChans = any(isnan(sigIn(:,:,:)),3);
    zeroChans = all(sigIn(:,:,:)==0, 3);
    badChannels = find(nanChans | zeroChans);
end
if ~exist('numChannels','var')
    numChannels = 1;
end
if ~exist('stdVal','var')
    stdVal = 2 ;
end
addpath(genpath([pwd,'/ToolOthers/SpikeNet']))
%% Initialisation

GammaBurstEvent.spike_sort_range = 6;
GammaBurstEvent.no_std = stdVal; % 2 std above the mean 1.5
GammaBurstEvent.peak_no_std = stdVal; % 4 std above the mean 3
GammaBurstEvent.burst_min_ms = min_time;

dt = 1/fsTemporal*1000;

% burst detection
steps = size(sigIn,2);

burst_min_steps = round(GammaBurstEvent.burst_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
GammaBurstEvent.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
GammaBurstEvent.Hz = zeros(1,numChannels);
GammaBurstEvent.is_burst = zeros(numChannels,size(sigIn,2));

iter_num = 100;
hil_mean_baseline_hist = zeros(numChannels,iter_num);
hil_std_baseline_hist = zeros(numChannels,iter_num);

%% find Gamma Burst
 for iChan = setdiff(1:size(sigIn,1), badChannels)      
     LFP_gamma_hilbert_tmp = sigIn(iChan,:) ;
     is_burst_tmp = zeros(size(LFP_gamma_hilbert_tmp));
     for r_iter = 1:iter_num
         hil_mean_baseline_hist(iChan,r_iter) =  mean(LFP_gamma_hilbert_tmp(~is_burst_tmp));
         hil_std_baseline_hist(iChan,r_iter)  =  std(LFP_gamma_hilbert_tmp(~is_burst_tmp));
         
         if r_iter > 1
             mean_diff = abs(hil_mean_baseline_hist(iChan,r_iter) - ...
                 hil_mean_baseline_hist(iChan,r_iter-1))/hil_mean_baseline_hist(iChan,r_iter-1);
             std_diff = abs(hil_std_baseline_hist(iChan,r_iter) - ...
                 hil_std_baseline_hist(iChan,r_iter-1))/hil_std_baseline_hist(iChan,r_iter-1);
             if mean_diff < 0.05 && std_diff < 0.05
                 GammaBurstEvent.is_burst(iChan,:) = is_burst_tmp;
                 break
             end
         end
         
         % thresholding
         is_burst_tmp = LFP_gamma_hilbert_tmp - hil_mean_baseline_hist(iChan,r_iter) ...
             > GammaBurstEvent.no_std*hil_std_baseline_hist(iChan,r_iter);
         
         % Discard transient data
         is_burst_tmp(1:transient_steps) = 0;
         
         % persistence_requirement
         is_burst_tmp = persistence_requirement( is_burst_tmp, burst_min_steps );
         
         
         % peak height requirement
         peak_min = hil_mean_baseline_hist(iChan,r_iter) + ...
             GammaBurstEvent.peak_no_std*hil_std_baseline_hist(iChan,r_iter);
         [ is_burst_tmp ] = peak_height_requirement( is_burst_tmp,...
             LFP_gamma_hilbert_tmp, peak_min );
         
     end
     cutoff = 1;
     [~, burst_du, flat_du, burst_start, ~] = seq_postprocess(GammaBurstEvent.is_burst(iChan,:), 1, cutoff);
     GammaBurstEvent.Hz(iChan) = length(burst_du)/(dt*steps*10^-3); % burst frequency in whole simulation
     GammaBurstEvent.burst_du_steps{iChan} = burst_du;
     GammaBurstEvent.flat_du_steps{iChan} = flat_du;
     GammaBurstEvent.burst_start_steps{iChan} = burst_start;
     GammaBurstEvent.mean_baseline{iChan} = hil_mean_baseline_hist ;
     % figure
     % %subplot(3,3,count);
     % count = count+1 ;
     % axisEnd = length(GammaBurstEvent.burst_start_steps{1}) ;
     % errorbar(1:axisEnd,GammaBurstEvent.burst_start_steps{1}(1:axisEnd)/fsTemporal,...
     %     zeros(1,axisEnd),GammaBurstEvent.burst_du_steps{1}(1:axisEnd)/fsTemporal,...
     %     zeros(1,axisEnd),zeros(1,axisEnd),'.')
     % ylabel('time (s)')
     % xlabel('burst count')
     % %xlim([0,180])
     % title(['Gamma Burst at electrode (',int2str(xPoint), ', ', int2str(yPoint),')' ])
 end

rmpath(genpath([pwd,'/ToolOthers/SpikeNet']))
