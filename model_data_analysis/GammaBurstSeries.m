function GammaBurstSeries
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
dt = 0.1;
fs = 1/(dt*1e-3); % sampling frequency (Hz)
% Butterworth filter
order = 4; % 4th order
lowFreq = 4; % gamma band
hiFreq = 7;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.

for i = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', i, num_files);
    fprintf('\t File name: %s\n', files{i});
    R = load(files{i});
%     R = GetBurst(R);
    %     disp('Finish getting burst...')
    %     duration = 0.1*[R.LFP.GammaBurstEvent.burst_du_steps{:}]; % ms
    %     interval = 0.1*[R.LFP.GammaBurstEvent.flat_du_steps{:}]; % ms
    %     [Dfre,Duration,~,~] = BurstDrifting(R);
    %     dfre = [Dfre{:}];
    %     subplot(2,2,1)
    %     histogram(duration,20)
    %     text(-0.18,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    %     xlabel('Duration(ms)')
    %     ylabel('Count')
    %     subplot(2,2,2)
    %     histogram(interval,20)
    %     text(-0.18,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    %     xlabel('Interval(ms)')
    %     ylabel('Count')
    %     disp('Finish duration and interval...')
    %     subplot(2,2,3)
    %     histogram(dfre,20)
    %     text(-0.18,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    %     xlabel('Central Frequency(Hz)')
    %     ylabel('Count')
    
%     IsBurst = sum(R.LFP.GammaBurstEvent.is_burst)/16;
%     IsBurst = sum(R.LFP.GammaBurstEvent.is_burst([1 11],:),1)/2;
%     BurstRate = smooth(mean(vec2mat(IsBurst,10),2),30);
%     plot(1:R.step_tot/10,BurstRate)
%     xlabel('Time(ms)')
%     ylabel('Burst Rate(Hz)')
LFP = R.LFP.LFP{1};
[no,~] = size(LFP);
LFP_theta_hilbert_abs = zeros(size(LFP));
gaus_width = 125; % ms
[ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );
LFP_theta = zeros(size(LFP)); 
for i = 1:no
    LFP_theta(i,:) = filter(b,a,LFP(i,:));
%     LFP_gamma_hilbert_abs(i,:) = abs(hilbert(LFP_gamma(i,:)));
%     % hilbert transformation & gaussian smoothing
    LFP_theta_hilbert_abs(i,:) = conv(abs(hilbert(LFP_theta(i,:))), Kernel,'same'); 
end
% plot((1:R.step_tot)*0.1,mean(R.LFP.LFP_gamma_hilbert_abs),(1:R.step_tot)*0.1,mean(LFP_theta_hilbert_abs))
% legend('Gamma','Theta')
% plot((1:R.step_tot)*0.1,R.LFP.LFP_gamma_hilbert_abs(1,:),(1:R.step_tot)*0.1,R.LFP.LFP_gamma_hilbert_abs(11,:),(1:R.step_tot)*0.1,R.LFP.LFP_gamma_hilbert_abs(3,:))
% legend('No.1','No.11','No.3')
% xlabel('Time(ms)')
% ylabel('Power')

max_lag = 100; % ms 
dt = 0.1;
[ac1,lags1] = xcorr(mean(R.LFP.LFP_gamma_hilbert_abs),mean(LFP_theta_hilbert_abs),round(max_lag/dt),'coeff');
[ac2,lags2] = xcorr(mean(R.LFP.LFP_gamma),mean(LFP_theta),round(max_lag/dt),'coeff');
plot(lags1*dt,ac1,'b',lags2*dt,ac2,'r')
legend('Power','LFP')
xlabel('Lags(ms)')
ylabel('Gamma-Theta xcorr')
    next = input('\t Next figure?');
    close all
end
% disp('Start CFC...')
% CFC
end