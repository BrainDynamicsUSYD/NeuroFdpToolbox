%% This function finds the maximum frequency in the wavelet spectrum
% also compares to the 1D burst detection
%
% author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
%%
flagBandstop = 1 ;
[dataFull,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
fsTemporal = Fs ;

%% find Gamma Burst
addpath([pwd,'/ToolOthers/SpikeNet'])
count = 1 ;
GammaBurstEvent.spike_sort_range = 6;
GammaBurstEvent.no_std = 2; % 2 std above the mean 1.5
GammaBurstEvent.peak_no_std = 2; % 4 std above the mean 3
GammaBurstEvent.burst_min_ms = 15;

dt = 1/fsTemporal*1000;

% burst detection
no = 1;  % number of electrodes
steps = size(dataFull,3);

burst_min_steps = round(GammaBurstEvent.burst_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
GammaBurstEvent.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
GammaBurstEvent.index1 = cell(100,no);
GammaBurstEvent.index2 = cell(100,no);
GammaBurstEvent.index3 = cell(100,no);
GammaBurstEvent.inside_rate = cell(100,no);
GammaBurstEvent.outside_rate = cell(100,no);
GammaBurstEvent.Hz = zeros(1,no);
GammaBurstEvent.is_burst = zeros(100,size(dataFull,3));

iter_num = 100;
hil_mean_baseline_hist = zeros(no,iter_num);
hil_std_baseline_hist = zeros(no,iter_num);

for yPoint = 6% 1:10
    for xPoint = 4% 1:10
        % xPoint = 5 ;
        % yPoint = 5 ;
        if ~isnan(dataFull(xPoint,yPoint,1))
            tempData = squeeze(dataFull(xPoint,yPoint,:)) ;
            N = 8;  % Filter order
            % discardTimeSecs = 5;
            % Filter LFPs
            disp('Filtering waveforms...') ; tic
            fLow = 25 ;
            fHigh = 35 ;
            % Design Butterworth band-pass filter
            h = fdesign.bandpass('N,F3dB1,F3dB2',N,fLow,fHigh,fsTemporal);
            Hd = design(h, 'butter');
            set(Hd, 'Arithmetic', 'double');
            
            SOS = Hd.sosMatrix;
            G = Hd.ScaleValues;
            
            % Filter signals forwards and backwards to avoid phase distortion
            bandpassLFPs = filtfilt(SOS,G,tempData);
            
            % i = 1;
            i = count ;
            LFP_gamma_hilbert_tmp = abs(hilbert(bandpassLFPs)) ;
            is_burst_tmp = zeros(size(LFP_gamma_hilbert_tmp));
            for r_iter = 1:iter_num
                hil_mean_baseline_hist(i,r_iter) =  mean(LFP_gamma_hilbert_tmp(~is_burst_tmp));
                hil_std_baseline_hist(i,r_iter)  =  std(LFP_gamma_hilbert_tmp(~is_burst_tmp));
                
                if r_iter > 1
                    mean_diff = abs(hil_mean_baseline_hist(i,r_iter) - ...
                        hil_mean_baseline_hist(i,r_iter-1))/hil_mean_baseline_hist(i,r_iter-1);
                    std_diff = abs(hil_std_baseline_hist(i,r_iter) - ...
                        hil_std_baseline_hist(i,r_iter-1))/hil_std_baseline_hist(i,r_iter-1);
                    if mean_diff < 0.05 && std_diff < 0.05
                        GammaBurstEvent.is_burst(i,:) = is_burst_tmp;
                        break
                    end
                end
                
                % thresholding
                is_burst_tmp = LFP_gamma_hilbert_tmp - hil_mean_baseline_hist(i,r_iter) ...
                    > GammaBurstEvent.no_std*hil_std_baseline_hist(i,r_iter);
                
                % Discard transient data
                is_burst_tmp(1:transient_steps) = 0;
                
                % persistence_requirement
                is_burst_tmp = persistence_requirement( is_burst_tmp, burst_min_steps );
                
                
                % peak height requirement
                peak_min = hil_mean_baseline_hist(i,r_iter) + ...
                    GammaBurstEvent.peak_no_std*hil_std_baseline_hist(i,r_iter);
                [ is_burst_tmp ] = peak_height_requirement( is_burst_tmp,...
                    LFP_gamma_hilbert_tmp, peak_min );
                
            end
            cutoff = 1;
            [~, burst_du, flat_du, burst_start, ~] = seq_postprocess(GammaBurstEvent.is_burst(i,:), 1, cutoff);
            GammaBurstEvent.Hz(i) = length(burst_du)/(dt*steps*10^-3); % burst frequency in whole simulation
            GammaBurstEvent.burst_du_steps{i} = burst_du;
            GammaBurstEvent.flat_du_steps{i} = flat_du;
            GammaBurstEvent.burst_start_steps{i} = burst_start;
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
            count = count + 1;
        end
    end
end


%% plot the power spectrum
% close all
for timeStep = 1
    addpath('ToolOthers/uimage')
    timeEpoch = 80 ;
    timeStart = timeEpoch*(timeStep-1)  ;
    timeEnd = timeEpoch*timeStep   ;
    %timeStart = 0 ;
    %timeEnd = 20 ;
    
    freqRange = 20:100 ;
    fc = centfrq('cmor1.5-1') ;
    scalerange = fc./(freqRange/fsTemporal) ;
    scales = scalerange(end):0.5:scalerange(1) ;
    pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
    
    tempData =  squeeze (dataFull(4,6,:)) ; %floor(timeStart*fsTemporal)+1 : ...
       %floor(timeEnd*fsTemporal))) ; 
%     tempData = squeeze(dataFull( 8,floor(timeStart*fsTemporal)+1 : ...
%         floor(timeEnd*fsTemporal) ))  ;
    
    wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
    tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
    
%     figure
%     wtNorm = zscore(abs(wt),[],2) ;
%     
%     imagesc(tempTimeAxis,pseudoFreq,(wtNorm(:,:)))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     figure
%     % uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
end

%%
freqCount = 1 ;
absWavelet = abs(wt) ;
for iBurst = 1:300%length(GammaBurstEvent.burst_start_steps{1})
    burstStart = GammaBurstEvent.burst_start_steps{1}(iBurst) ;
    burstDu = GammaBurstEvent.burst_start_steps{1}(iBurst) ;
    
    for iTime = burstStart:burstStart + burstDu
        [~, maxFreqIdx] = max(absWavelet(:,iTime)) ;
        maxFreq(freqCount) = pseudoFreq(maxFreqIdx) ;
        freqCount = freqCount+1 ;
    end
    
end
%%

hist(maxFreq(1:100000),200)
title('Maximum frequency distribution for 100s for single electrode')


%%
maxFreqTemp = maxFreq(maxFreq>=30) ; 
[x,n] = hist(maxFreqTemp,80) ;
loglog(n,x/sum(x),'o')
hist(maxFreqTemp,40) ;
title('max frequency within bursts (100 bursts)')
xlabel('frequency (Hz)')
ylabel('probability')

xmin = 31 ;
sigIn = maxFreq ;
% sigIn = sigIn(sigIn>xmin) ;
%%
pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
[lambdaHat1,lambdaCI] = mle(maxFreqTemp, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
L1 = sum(log(pf_power(maxFreqTemp,lambdaHat1 ))) ;

%%
PARMHAT = lognfit(maxFreqTemp) ;
[n,x] = hist(maxFreqTemp,40) ;
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
title('Maximum frequency distribution')
xlabel('freq')
ylabel('probability')


%% log-normal Burst interval (single electrode)

PARMHAT = lognfit(flat_du) ;
[n,x] = hist(flat_du,40) ;
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
%%
hist(flat_du,40)
title('Burst interval (single electrode)')
xlabel('time (ms)')
ylabel('probability')

%% power law
figure;
loglog(x,n/sum(n),'o')

%% log-normal Burst interval (3D)
PARMHAT = lognfit(centInterval(centInterval>0)) ;
[n,x] = hist(centInterval(centInterval>0),80) ;
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
title('burst interval distribution')
xlabel('time(ms)')
ylabel('probability')

%%
hist(GammaBurstEvent.burst_du_steps{1},30)
title('burst duration distribution')
xlabel('time(ms)')
ylabel('count')

%% log-normal Burst interval (3D)
PARMHAT = lognfit(GammaBurstEvent.burst_du_steps{1}) ;
[n,x] = hist(GammaBurstEvent.burst_du_steps{1},40) ;
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
title('burst duration distribution')
xlabel('time(ms)')
ylabel('probability')


%% log-normal width

a = []; 
for itemp = 1:3272;
a=[a;width{itemp}];
end
hist(a,18)
xlabel('electrode')
ylabel('count')
title('distribution of width')
%%

PARMHAT = lognfit(a) ;
[n,x] = hist(a,20) ;
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
title('burst duration distribution')
xlabel('time(ms)')
ylabel('probability')

